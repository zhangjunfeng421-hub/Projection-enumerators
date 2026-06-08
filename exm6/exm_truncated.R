source("PSE.R")
source("PBE.R")

# Construct OSOA(2m, m-1, 4, 3-) by Hadamard matrix according to Zhou and Tang (2019)
OSOA <- function(m) {
  if (m %% 4 != 0) {
    return(NULL)
  } else {
    H <- hadamard(m)
    H_1 <- cbind(H, H)
    H_2 <- cbind(H, -H)
    HH <- rbind(H_1, H_2)
    A <- HH[, (m + 2):(2 * m)]
    A <- (A + 1) / 2
    B <- HH[, 2:m]
    B <- (B + 1) / 2
    SOA <- 2 * A + B
    return(SOA)
  }
}

A <- OSOA(96)
time_record1 <- list()
dims_truncated <- c(2,5,10,20)
for (i in 1:4) {
  time_start <- Sys.time()
  TPSFPI(A,2,2,dims_truncated[i])
  time_end <- Sys.time()
  time_record1[[i]] <- as.numeric(time_end - time_start, units = "secs")
  print(time_record1[[i]])
}
time_start <- Sys.time()
pattern_a <- PSFPI(A,2,2)
time_end <- Sys.time()
time_record1[[5]] <- as.numeric(time_end - time_start, units = "secs")
print(time_record1[[5]])
time_record1
save(time_record1, file = "time_TPSE.Rdata")

time_record2 <- list()

# Generator matrix
element_OA <- function(s,p){
  G <- matrix(0:(s - 1), ncol = 1)
  l <- matrix(1, nrow = s, ncol = 1)
  if(p == 1){
    return(G)
  }
  
  for (i in 1:(p - 1)) {
    a <- G[, 1] %x% l
    G <- cbind(a, l %x% G)
  }
  return(G)
}


## A function to generate a saturated orthogonal table. 
OA_Saturated <- function(s, p){
  G <- element_OA(s,p)
  OA <- matrix(0, nrow = s^p, ncol = (s^p - 1) / (s - 1))
  t <- 0
  
  for (i in 1:p) {
    if(i == 1){
      OA[, 1] <- G[, 1]
    }else{
      P <- t(element_OA(s, (i - 1)))
      OA[, (t + 1): (t + s^(i - 1))] <- ((G[, 1:(i - 1)] %*% P + G[, i]) %% s)
    }
    t <- t + s^(i - 1)
  }
  return(OA)
}
A <- OA_Saturated(3,5)
dims_truncated <- c(2,5,10,20)
for (i in 1:4) {
  time_start <- Sys.time()
  TPWPI(A,3,dims_truncated[i])
  time_end <- Sys.time()
  time_record2[[i]] <- as.numeric(time_end - time_start, units = "secs")
  print(time_record2[[i]])
}
time_start <- Sys.time()
pattern_b <- PWPI(A,3)
time_end <- Sys.time()
time_record2[[5]] <- as.numeric(time_end - time_start, units = "secs")
time_record2
save(time_record2, file = "time_TPBE.Rdata")
