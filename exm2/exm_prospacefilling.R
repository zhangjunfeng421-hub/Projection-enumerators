library(foreach)
library(doParallel)
source("PSE.R")
source("stratification_pattern.R")

# Generate a design according to Shi and Tang (2020), calculate its PSFP and 
# record the calculation time 
A <- SOA_tang(7)
time_record <- list(NULL)
maxm <- 5
v <- seq(1,31)

for (m in 1:9) {
  print(paste("m is", m))
  

  num_cores <- 12
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  clusterExport(cl, c("A", "v", "m", "maxm", "D_SP"))
  clusterEvalQ(cl, {
    source("PSE.R")
    source("stratification_pattern.R")
  })
  
  time_results <- foreach(i = 1:maxm, 
                          .combine = 'c') %dopar% {
                            colname <- sample(v, size = m, replace = FALSE)
                            D <- A[, colname, drop = FALSE]  
                            
                            time_start <- Sys.time()
                            P <- D_SP(D, 3, 3 * m)
                            time_end <- Sys.time()
                            time_execution <- time_end - time_start
                            
                            print(paste("m =", m, "i =", i, "time =", time_execution))
                            
                            return(as.numeric(time_execution))
                          }
  
  stopCluster(cl)
  
  time_record[[m]] <- mean(time_results)
  print(paste("m =", m, "average time =", time_record[[m]]))
}

print(time_record)
save(time_record, file = "exm_prospacefilling.Rdata")
