library(Matrix)
library(pracma)
library(gtools)

# Define an orthogonal polynomial function; we suggest N <= 88;
#for the enumerator appropriating more levels, 
#one can use a more stable algorithm to generate the similarity contrast(the key is to stable polynomial_contrast())
polynomial_contrast <- function(N) {
  if(N == 1){
    return(matrix(1,1,1))
  }
  x <- 0:(N-1)
  
  phi <- matrix(0, nrow = N, ncol = N)
  
  phi[,1] <- rep(1, N)
  phi[,1] <- phi[,1] * sqrt(N / sum(phi[,1]^2))
  
  for (k in 2:N) {
    v <- x^(k-1)
    
    for (j in 1:(k-1)) {
      proj_coef <- sum(v * phi[,j]) / N
      v <- v - proj_coef * phi[,j]
    }
    norm2 <- sum(v^2)
    phi[,k] <- v * sqrt(N / norm2)
  }
  
  return(phi)
}

# To define a function to calculate the PBE, we need to calculate the similarity contrast
Contrast_similarity_array <- function(N){
 A <- polynomial_contrast(N)
 result <- array(0,dim = c(N,N,N))
 for (k in 1:N) {
   result[, , k] <- outer(A[, k], A[, k])
 }
 return(result)
}

# Construct the table to record the R_varphi between any two levels
build_R_varphi_table <- function(y, contrast_similarity) {
  N <- dim(contrast_similarity)[1]
  v <- y^(0:(N-1))

  mat <- matrix(contrast_similarity, nrow = N*N, ncol = N)
  
  result_vec <- mat %*% v
  
  result <- matrix(result_vec, nrow = N, ncol = N)
  
  return(result)
}

# Define a function to calculate the PBE
PBE <- function(x, y, D, s) {
  Contrast_similarity <- Contrast_similarity_array(s)
  R_table <- build_R_varphi_table(y,Contrast_similarity)
  W <- 1 + x * (R_table - 1)
  L <- log(W)                         
  D_idx <- D + 1
  n <- nrow(D)
  m <- ncol(D)
  F <- matrix(0, n, n)
  for (k in 1:m) {
    col_vals <- D_idx[, k]
    F <- F + L[col_vals, col_vals]
  }
  PBE_val <- mean(exp(F))
  return(PBE_val)
}

# Calculate PBE given the similarity contrast
PBE_similarity <- function(x, y, D, Contrast_similarity) {
  R_table <- build_R_varphi_table(y,Contrast_similarity)
  W <- 1 + x * (R_table - 1)
  L <- log(W)                         
  D_idx <- D + 1
  n <- nrow(D)
  m <- ncol(D)
  F <- matrix(0, n, n)
  for (k in 1:m) {
    col_vals <- D_idx[, k]
    F <- F + L[col_vals, col_vals]
  }
  PBE_val <- mean(exp(F))
  return(PBE_val)
}

# Using PBE to Inversely Solve \beta-PWP
PWPI <- function(D, s) {
  n <- nrow(D)
  m <- ncol(D)
  Similarity <- Contrast_similarity_array(s)
  Ep <- matrix(0, nrow = (m * (s-1)), ncol = m)
  
  for (k in 1:(m * (s-1))) {
    y <- exp((1i * k * 2 * pi) / (m * (s-1)))
    for (j in 1:m) {
      x <- exp((1i * j * 2 * pi) / m)
      Ep[k, j] <- PBE_similarity(x, y, D, Similarity)
    }
  }
  Ep <- Ep - 1
  Ep_fft <- fft(Ep)
  B_matrix <- Mod(Ep_fft) / (m^2 * (s - 1))
  B_matrix <- cbind(B_matrix[,2:m],B_matrix[,1])
  B_matrix <- rbind(B_matrix[2:(m*(s-1)),],B_matrix[1,])

  indices <- expand.grid(i = 1:(m * (s-1)), j = 1:m)
  indices <- indices[indices$j >= ceiling(indices$i / (s-1)) &
                       indices$j <= indices$i, ]
  indices <- indices[order(indices$i, indices$j), ]
  B_vector <- B_matrix[cbind(indices$i, indices$j)]
  B_vector <- round(B_vector, 4)
  return(B_vector)
}
# A function to calculate truncated PBE
TPBE  <- function(x, y, D,s,t) {
  Contrast_similarity <- Contrast_similarity_array(s)
  R_table <- build_R_varphi_table(y,Contrast_similarity)
  R_table <- x * (R_table - 1) 
  n <- nrow(D)
  m <- ncol(D)
  
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
# Using similarity contrast to calculate R_table for y by symmetrical elementary polynomial 
elementary_result_TBE <- function(y, D, Contrast_similarity,t) {
  R_table <- build_R_varphi_table(y,Contrast_similarity)
  R_table <- R_table - 1
  n <- nrow(D)
  m <- ncol(D)
  
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
# Using TPBE to Inversely Calculate $\beta$-PWP restricted dims-dimension
TPWPI <- function(D, s, dims) {
  n <- nrow(D)
  m <- ncol(D)
  Similarity <- Contrast_similarity_array(s)
  Ep <- matrix(0, nrow = (dims * (s-1)), ncol = dims)
  
  for (k in 1:(dims * (s-1))) {
    y <- exp((1i * k * 2 * pi) / (dims * (s-1)))
    S_y <- elementary_result_TBE(y, D, Similarity, dims)
    for (j in 1:dims) {
      x <- exp((1i * j * 2 * pi) / dims)
      X <- array(rep(x^(0:dims), each = n * n), dim = c(n, n, dims + 1))
      S_xy <- X * S_y
      S_xy <- rowSums(S_xy, dims = 2)
      Ep[k, j] <- mean(S_xy) 
    }
  }
  Ep <- Ep - 1
  B_matrix <- fft(Ep) / (dims * (dims * (s - 1)))
  B_matrix <- cbind(B_matrix[,2:dims],B_matrix[,1])
  B_matrix <- rbind(B_matrix[2:(dims*(s - 1)),],B_matrix[1,])
  
  indices <- expand.grid(i = 1:(dims * (s - 1)), j = 1:dims)
  indices <- indices[indices$j >= ceiling(indices$i / (s - 1)) & 
                       indices$j <= indices$i, ]
  indices <- indices[order(indices$i, indices$j), ]
  B_vector <- Mod(B_matrix[cbind(indices$i, indices$j)])
  B_vector <- round(B_vector, 4)
  
  return(B_vector)
}