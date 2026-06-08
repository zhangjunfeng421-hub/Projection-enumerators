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

# projection stratification pattern of D_1
PSP <- PSFPI(D_1,2,3)
save(PSP, file = "PSP.Rdata")

# $\beta$-projection wordlength pattern of D_2
PWP_beta <- PWPI(D_2,8)
save(PWP_beta, file = "PWP.Rdata")
