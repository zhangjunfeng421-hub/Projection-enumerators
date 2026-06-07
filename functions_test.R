# ============================================================
# Define the real functions and the ranges of the explanatory variables
# ============================================================
# 1. Borehole function
Borehole <- function(x) {
  
  x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
  x5 <- x[5]; x6 <- x[6]; x7 <- x[7]; x8 <- x[8]
  
  numerator <- 2 * pi * x3 * (x4 - x6)
  denominator <- log(x2 / x1) * (1 + (2 * pi * x3) / (log(x2 / x1) * x1^2 * x8) + x3 / x5)
  
  y <- numerator / denominator
  return(y)
}

Ranges_Borehole <- matrix(c(
  0.05, 0.15,      # x1
  100, 50000,      # x2
  63070, 115600,   # x3
  990, 1110,       # x4
  63.1, 116,       # x5
  700, 820,        # x6
  1120, 1680,      # x7
  9855, 12045      # x8
), nrow = 8, ncol = 2, byrow = T)

# 2. Wing weight functions
Wing <- function(x) {
  S_w <- x[1]      # wing area (ft^2)
  W_fw <- x[2]     # weight of fuel in the wing (lb)
  A <- x[3]        # aspect ratio
  Lambda <- x[4]   # quarter-chord sweep (degrees)
  q <- x[5]        # dynamic pressure at cruise (lb/ft^2)
  lambda <- x[6]   # taper ratio
  t_c <- x[7]      # aerofoil thickness to chord ratio
  N_z <- x[8]      # ultimate load factor
  W_dg <- x[9]     # flight design gross weight (lb)
  W_p <- x[10]     # paint weight (lb/ft^2)
  
  term1 <- 0.036 * (S_w^0.758) * (W_fw^0.0035) * 
    ((A / (cos(Lambda * pi / 180)^2))^0.6) * 
    (q^0.006) * (lambda^0.04) * 
    ((100 * t_c / cos(Lambda * pi / 180))^(-0.3)) * 
    ((N_z * W_dg)^0.49)
  
  term2 <- S_w * W_p
  
  y <- term1 + term2
  return(y)
}

Ranges_Wing <- matrix(c(
  150, 200,      # x1
  220, 300,      # x2
  6, 10,   # x3
  -10, 10,       # x4
  16, 45,       # x5
  0.5, 1,        # x6
  0.08, 0.18,      # x7
  2.5, 6,     # x8
  1700, 2500,     # x9
  0.025, 0.08     # x10
), nrow = 10, ncol = 2, byrow = T)

# 3. OTL circuit function (Ben-Ari & Steinberg, 2007)
OTL <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  x5 <- x[5]
  x6 <- x[6]
  
  term1 <- (12 * x2 / (x1 + x2) + 0.74) * x6 * (x5 + 9)
  term2 <- (11.35 * x3) / (x6 * (x5 + 9) + x3)
  term3 <- (0.74 * x3 * x6 * (x5 + 9)) / ((x6 * (x5 + 9) + x3) * x4)
  
  y <- term1 + term2 + term3
  return(y)
}

Ranges_OTL <-  matrix(c(
  50, 150,      # x1
  25, 70,      # x2
  0.5, 3,      # x3
  1.2, 2.5,      # x4
  0.25, 1.2,      # x5
  50, 300       # x6
), nrow = 6, ncol = 2, byrow = T)

# Piston simulation function
Piston <- function(x) {
  
  M <- x[1]
  k <- x[2]
  S <- x[3]
  P0 <- x[4]
  V0 <- x[5]
  T0 <- x[6]
  T <- x[7]
  
  A <- P0 * S + 19.62 * M - (k * V0) / S
  
  V <- (S / (2 * k)) * (sqrt(A^2 + 4 * k * (P0 * V0 / T0) * T) - A)
  
  cycle_time <- 2 * pi * sqrt(M / (k + S^2 * (P0 * V0 / T0) * (T / V^2)))
  
  return(cycle_time)
}

Ranges_Piston <-  matrix(c(
  30, 60,     # x1
  0.005, 0.020, # x2
  0.002, 0.010,      # x3
  1000, 5000,      # x4
  90000, 110000,      # x5
  290, 296,       # x6
  340, 360        # x7
), nrow = 7, ncol = 2, byrow = T)
# 5. MOON ET AL. (2012) FUNCTION
#set.seed(2026)
#n_samples <- 1000000  

#X_bh <- matrix(runif(n_samples * 8), n_samples, 8)
#for (j in 1:8) {
#  X_bh[, j] <- Ranges_Borehole[j, 1] + X_bh[, j] * 
#    (Ranges_Borehole[j, 2] - Ranges_Borehole[j, 1])
#}
#Y_bh <- apply(X_bh, 1, Borehole)

#X_wing <- matrix(runif(n_samples * 10), n_samples, 10)
#for (j in 1:10) {
#  X_wing[, j] <- Ranges_Wing[j, 1] + X_wing[, j] * 
#    (Ranges_Wing[j, 2] - Ranges_Wing[j, 1])
#}
#Y_wing <- apply(X_wing, 1, Wing)

#X_otl <- matrix(runif(n_samples * 6), n_samples, 6)
#for (j in 1:6) {
#  X_otl[, j] <- Ranges_OTL[j, 1] + X_otl[, j] * 
#    (Ranges_OTL[j, 2] - Ranges_OTL[j, 1])
#}
#Y_otl <- apply(X_otl, 1, OTL)

#X_piston <- matrix(runif(n_samples * 7), n_samples, 7)

#for (j in 1:7) {
#  X_piston[, j] <- Ranges_Piston[j, 1] + X_piston[, j] * 
#    (Ranges_Piston[j, 2] - Ranges_Piston[j, 1])
#}

#Y_piston <- apply(X_piston, 1, Piston)

#maxy1 <- max(Y_bh); maxy2 <- max(Y_wing); maxy3 <- max(Y_otl); maxy4 <- max(Y_piston)
#miny1 <- min(Y_bh); miny2 <- min(Y_wing); miny3 <- min(Y_otl); miny4 <- min(Y_piston)

#Moon <- function(x){
#  x1 <- x[1:8]; x2 <- x[9:18]; x3 <- x[19:24]; x4 <- x[25:31];
#  y1 <- Borehole(x1); y2 <- Wing(x2); y3 <- OTL(x3); y4 <- Piston(x4)
#  y <- (y1 - miny1)/(maxy1- miny1) + (y2 - miny2)/(maxy2- miny2) + (y3 - miny3)/(maxy3- miny3) + (y4 - miny4)/(maxy4- miny4)
 # return(y)
#}
#Ranges_Moon <- rbind(rbind(rbind(Ranges_Borehole,Ranges_Wing), Ranges_OLT), Ranges_Piston)