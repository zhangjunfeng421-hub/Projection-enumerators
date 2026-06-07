library(Matrix)
library(pracma)
library(gtools)

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

# Define a function to calculate the PSE given distance matrix D of a design 
R_chi <- function(y, rho, s, p) {
  
  if (rho == 0) {
    b <- ((1 - y) * (1 - (s * y)^(p + 1)) / (1 - s * y)) + (s^p) * (y^(p + 1)) 
  } else {
    b <- ( (1 - y) * (1 - (s * y)^(p + 1 - rho)) / (1 - s * y)) 
  }
  
  return(b)
}

# Build NRT-distance matrix for 0,...,s^p-1
NRT_matrix <- function(s, p){
  S <- s^p
  all_vals <- 0:(S - 1)
  coord_all <- cbind(coordinate(all_vals, s, p),0)
  distance_matrix <- matrix(0, S, S)
  for (i in 1:S) {
    diff <- ((matrix(rep(coord_all[i,], s^p), nrow = s^p ,byrow = T) - coord_all) != 0)
    first_nonzero <- max.col(diff, ties.method = "first")
    first_nonzero[ rowSums(diff) == 0 ] <- p + 1
    distance_matrix[i, ] <- p + 1 - first_nonzero
  }
  
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

# Define a function to calculate the PSE
PSE <- function(x, y, D, s, p) {
  S <- s^p
  Rho <- NRT_matrix(s,p)
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
PSFPI <- function(D, s, p) {
  n <- nrow(D)
  m <- ncol(D)
  distance_matrix <- NRT_matrix(s,p)
  Ep <- matrix(0, nrow = (m * p), ncol = m)

  for (k in 1:(m * p)) {
    for (j in 1:m) {
      y <- exp((1i * k * 2 * pi) / (m * p))
      x <- exp((1i * j * 2 * pi) / m)
      Ep[k, j] <- PSE_distance(x, y, distance_matrix, D, s, p)
    }
  }
  Ep <- Ep - 1

  Ep_fft <- fft(Ep)
  P_matrix <- Mod(Ep_fft) / (m^2 * p)

  P_matrix <- cbind(P_matrix[,2:m],P_matrix[,1])
  P_matrix <- rbind(P_matrix[2:(m*p),],P_matrix[1,])
  indices <- expand.grid(i = 1:(m * p), j = 1:m)
  indices <- indices[indices$j >= ceiling(indices$i / p) &
                       indices$j <= indices$i, ]
  indices <- indices[order(indices$i, indices$j), ]
  P_vector <- P_matrix[cbind(indices$i, indices$j)]
  P_vector <- round(P_vector, 4)
  return(P_vector)
}

# A function to calculate truncated PSE
TPSE <- function(x, y, D, s, p, t) {
  n <- nrow(D)
  m <- ncol(D)
  S <- s^p
  
  Rho <- NRT_matrix(s, p)                   
  R_table <- build_R_table(s, p, y, rho = Rho) 
  R_table <- x * (R_table - 1)           
  D_idx <- D + 1 
  
  C <- array(0, dim = c(n, n, t + 1))
  C[,,1] <- 1                             
  
  for (k in 1:m) {
    col_vals <- D_idx[, k]                    
    R_k <- R_table[col_vals, col_vals]
    R_k <- array(R_k, dim = c(n, n, t))
    C[,,2:(t+1)] <- C[,,2:(t+1)] + R_k * C[,,1:t]
  }
  
  S_mat <- rowSums(C, dims = 2)       
  result <- mean(S_mat)
  return(result)
}

# Using distance matrix to calculate R_table for y by symmetrical elementary polynomial 
elementary_result  <- function(y, D, Rho,s, p, t) {
  n <- nrow(D)
  m <- ncol(D)
  S <- s^p

  R_table <- build_R_table(s, p, y, rho = Rho) 
  R_table <- R_table - 1           
  
  D_idx <- D + 1 
  
  C <- array(0, dim = c(n, n, t + 1))
  C[,,1] <- 1

  for (k in 1:m) {
    col_vals <- D_idx[, k]
    R_k <- R_table[col_vals, col_vals]
    R_k <- array(R_k, dim = c(n, n, t))
    C[,,2:(t+1)] <- C[,,2:(t+1)] + R_k * C[,,1:t]
  }

  # S_mat <- rowSums(C, dims = 2)
  # result <- mean(S_mat)
  return(C)
}

# Using TPSE to Inversely Calculate PSFP restricted dims-dimension
TPSFPI <- function(D, s, p, dims) {
  n <- nrow(D)
  m <- ncol(D)
  distance_matrix <- NRT_matrix(s,p)
  Ep <- matrix(0, nrow = (dims * p), ncol = dims)
  
  for (k in 1:(dims * p)) {
    y <- exp((1i * k * 2 * pi) / (dims * p))
    S_y <- elementary_result(y, D, distance_matrix, s, p, dims)
    for (j in 1:dims) {
      x <- exp((1i * j * 2 * pi) / dims)
      X <- array(rep(x^(0:dims), each = n * n), dim = c(n, n, dims + 1))
      S_xy <- X * S_y
      S_xy <- rowSums(S_xy, dims = 2)
      Ep[k, j] <- mean(S_xy) 
    }
  }
  Ep <- Ep - 1
  Ep_fft <- fft(Ep)
  P_matrix <- Mod(Ep_fft) / (dims^2 * p)
  
  P_matrix <- cbind(P_matrix[,2:dims],P_matrix[,1])
  P_matrix <- rbind(P_matrix[2:(dims * p),],P_matrix[1,])
  
  indices <- expand.grid(i = 1:(dims * p), j = 1: dims)
  indices <- indices[indices$j >= ceiling(indices$i / p) &
                       indices$j <= indices$i, ]
  indices <- indices[order(indices$i, indices$j), ]
  P_vector <- P_matrix[cbind(indices$i, indices$j)]
  P_vector <- round(P_vector, 4)
  
  return(P_vector)
}