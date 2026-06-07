#=====================================================
# Select a good design under $L_2$-distance criterion
source("GSOA_opt.R")
euclidean_distance <- function(D) {
  n <- nrow(D)
  min_dist <- Inf
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      dist <- sum((D[i,] - D[j,])^2)
      if(dist < min_dist) min_dist <- dist
    }
  }
  return(min_dist)
}
# Define an element, which is represented by a polyomial of p-1 degrees, to shift the GSOA_opt.
# So we use a p-dimension vector to represent it, and use a s^p \times p matrix to list all possible vectors
s = 3
p = 3
element_shift <- as.matrix(expand.grid(rep(list(0:(s-1)), p)))

Design_opt <- GSOA_opt(c(1,0,2,1), s = 3 , K = 1)
distance_before <- euclidean_distance(Design_opt)
shift <- c(0,0,0)
for (i in 2:(s^p)) {
  Design_comp <- GSOA_opt(c(1,0,2,1), s = 3, K = 1, add = element_shift[i,])
  distance_after <- euclidean_distance(Design_comp)
  if(distance_after > distance_before){
    distance_before <- distance_after
    Design_opt <- Design_comp
    shift <- element_shift[i,]
  }
}
print(shift)
