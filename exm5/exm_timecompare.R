library(Matrix)
library(pracma)
library(gtools)
library(parallel)
library(foreach)
library(doParallel)

# Compute the coefficients under s-base
coordinate <- function(u, s, p) {
  if (length(u) == 1) {
    powers <- s^((p-1):0)
    return((u %/% powers) %% s)
  } else {
    powers <- matrix(s^((p-1):0), nrow = length(u), ncol = p, byrow = TRUE)
    return((u %/% powers) %% s)
  }
}


# Define a function to compute the NRT-distance 
NRT <- function(coord_u, coord_v, s, p) {
  coord_u <- cbind(coord_u,0)
  coord_v <- cbind(coord_v,0)
  coord_different <- ((coord_u- coord_v) != 0)
  a <- max.col(coord_different,ties.method = "first")
  a[rowSums(coord_different) == 0] <- p + 1
  a <- p + 1 - a
  return(a)
}

# Build NRT-distance matrix for 0,...,s^p-1
NRT_matrix <- function(s, p, n_cores = detectCores() - 1 ){
  S <- s^p
  all_vals <- 0:(S - 1)
  coord_all <- cbind(coordinate(all_vals, s, p),0)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl) 
  
  distance_matrix <- foreach(i = 1:S, .combine = rbind, .export = c("coord_all", "s", "p")) %dopar% {
    diff <- ((matrix(rep(coord_all[i,], s^p), nrow = s^p ,byrow = T) - coord_all) != 0)
    first_nonzero <- max.col(diff, ties.method = "first")
    first_nonzero[ rowSums(diff) == 0 ] <- p + 1
    distance <- p + 1 - first_nonzero
  }
  
  stopCluster(cl)
  return(distance_matrix)
} 

# Record the similarity coNRTast in a table
build_R_table <- function(s, p, y, rho) {
  S <- s^p
  R_table <- matrix(0, S, S)
  
  idx0 <- (rho == 0)
  R_table[idx0] <- ((1 - y)*(1 - (s * y)^(p + 1))/(1 - s * y)) + (s^p) * (y^(p + 1))
  
  idx_pos <- !idx0
  R_table[idx_pos] <- (1-y)*(1-(s*y)^(p+1 - rho[idx_pos]))/(1-s*y)
  return(R_table)
}

# Calculate PSE based the NRT-distance matrix 
PSE_distance <- function(x, y, Rho, D, s, p) {
  S <- s^p
  R_table <- build_R_table(s, p, y, rho = Rho)  
  W <- 1 + x * (R_table - 1)
  D_idx <- D + 1
  n <- nrow(D)
  m <- ncol(D)
  F <- matrix(1, n, n)
  for (k in 1:m) {
    col_vals <- D_idx[, k]
    F <- F * W[col_vals, col_vals]
  }
  PSE_val <- mean(F)
  return(PSE_val)
}

#Using PSE to Inversely Solve PSFP
PSFPI <- function(D, s, p, ncores = detectCores() - 1) {
  n <- nrow(D)
  m <- ncol(D)
  distance_matrix <- NRT_matrix(s,p, ncores)
  
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl) 
  
  
  Ep <- foreach(k = 1:(m * p), .combine = rbind, .export = c("distance_matrix","D","s","p","m","build_R_table","PSE_distance")) %dopar%{
    y <- exp((1i * k * 2 * pi) / (m * p))
    E_vector <- rep(0,m)
    for (j in 1:m) {
      x <- exp((1i * j * 2 * pi) / m)
      E_vector[j] <- PSE_distance(x, y, distance_matrix, D, s, p)
    }
    E_vector
  }
  stopCluster(cl)
  Ep <- Ep - 1
  
  Ep_fft <- fft(Ep)
  P_matrix <- Mod(Ep_fft) / (m^2 * p)
  
  if (m > 1) {
    P_matrix <- cbind(P_matrix[, 2:m], P_matrix[, 1])
    P_matrix <- rbind(P_matrix[2:(m*p),], P_matrix[1,])
  } else {
    P_matrix <- rbind(P_matrix[2:(m*p),], P_matrix[1,])
  }
  if(m == 1){return(P_matrix)}
  indices <- expand.grid(i = 1:(m * p), j = 1:m)
  indices <- indices[indices$j >= ceiling(indices$i / p) &
                       indices$j <= indices$i, ]
  indices <- indices[order(indices$i, indices$j), ]
  P_vector <- P_matrix[cbind(indices$i, indices$j)]
  P_vector <- round(P_vector, 4)
  return(P_vector)
}

all_column_products <- function(mat) {
  n <- ncol(mat)
  k <- nrow(mat)
  result <- list()
  
  for (i in 1:n) {
    combinations <- combinations(n, i)
    for (j in 1:nrow(combinations)) {
      cols <- combinations[j, ]
      prod_mat <- apply(mat[, cols, drop = FALSE], 1, prod)
      result <- c(result, list(prod_mat))
    }
  }
  
  Result <- matrix(0, nrow = k, ncol = length(result))
  for (i in 1:length(result)) {
    Result[, i] <- result[[i]]
  }
  return(Result)
}

SOA_tang <- function(k){
  E <- matrix(c(1, -1), nrow = 2, ncol = 1, byrow = T)
  if(k <= 2){
    return(NULL)
  } else {
    for (j in 2:k){
      e_1 <- c(rep(1, each = 2^(j - 1)), rep(-1, each = 2^(j - 1)))
      e_1 <- as.matrix(e_1)
      E <- rbind(E, E)
      E <- cbind(E, e_1)
    }
    P <- as.matrix(E[, 3:ncol(E)])
    P <- all_column_products(P)
    e_1 <- E[, 1]
    e_2 <- E[, 2]
    e_3 <- E[, 3]
    A <- e_1 * P
    B <- e_2 * P
    p <- ncol(P)
    if(p == 1){
      C <- P
    } else {
      e <- as.array(e_1 * e_2 * e_3)
      C <- matrix(rep(e, each = p), ncol = p, byrow = T)
    }
    S <- 2 * A + B + C/2 + 7/2
    return(S)
  }
}


A <- SOA_tang(7)

time_record <- list(NULL)
max_repeats <- 100         
v <- seq(1, 25)           
num_cores <- 12           

for (m in 1:25) {
  time_results <- numeric(max_repeats)
  for (i in 1:max_repeats) {
    colname <- sample(v, size = m, replace = FALSE)
    D <- A[, colname, drop = FALSE]
    
    time_start <- Sys.time()
    P <- PSFPI(D, s = 2, p = 3, ncores = num_cores)
    time_end <- Sys.time()
    
    time_results[i] <- as.numeric(difftime(time_end, time_start, units = "secs"))
    cat(paste("  m =", m, "repetition", i, "/", max_repeats, 
              "time =", round(time_results[i], 4), "秒\n"))
  }
  time_record[[m]] <- mean(time_results)
}

for (m in 1:25) {
  cat(paste("m =", sprintf("%2d", m), ":", 
            sprintf("%8.4f", time_record[[m]]), "秒\n"))
}

save(time_record, file = "time_theorem_parallel.Rdata")
