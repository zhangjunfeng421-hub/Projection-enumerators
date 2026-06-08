source("PSE.R")
source("PBE.R")
D_1 <- c(1,2,3,0,7,4,5,6,
         1,3,7,5,2,0,4,6,
         0,6,2,4,3,5,1,7)
D_1 <- matrix(D_1,nrow = 3,byrow = T)
D_1 <- t(D_1)

D_2 <- matrix(c(0.5,1.5,2.5,3.5,-3.5,-2.5,-1.5,-0.5,
                #-1.5,0.5,-3.5,2.5,-2.5,3.5,-0.5,1.5,
                3.5,2.5,-1.5,-0.5,0.5,1.5,-2.5,-3.5,
                2.5,-3.5,-0.5,1.5,-1.5,0.5,3.5,-2.5),nrow = 3,byrow = TRUE)
D_2 <- t(D_2)+3.5

Y <- c(1,0.1,0.01,0.001)
# projection stratification pattern of D_1
PSP_1 <- PSFPI(D_1,2,3)
PSP_2 <- PSFPI(D_2,2,3)
save(PSP_1, file = "PSP1.Rdata")
save(PSP_2, file = "PSP2.Rdata")

# $\beta$-projection wordlength pattern of D_2
PWP_1 <- PWPI(D_1,8)
save(PWP_1, file = "PWP1.Rdata")
PWP_2 <- PWPI(D_2,8)
save(PWP_2, file = "PWP2.Rdata")

PSE_value <- rep(0,4)
PBE_value <- rep(0,4)
Similarity <- Contrast_similarity_array(8)
for (i in 1:4) {
  y <- Y[i]
  PSE_value[i] <- PSE(y^(1/3), y, D_1, 2, 3)-PSE(y^(1/3), y, D_2, 2, 3)
  PBE_value[i] <- PBE(y^(1/3), y, D_1, Similarity)-PBE(y^(1/3), y, D_2, Similarity)
}
save(PSE_value, file = "PSE_difference.Rdata")
save(PBE_value, file = "PBE_difference.Rdata")
