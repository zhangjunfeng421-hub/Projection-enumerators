source("GSOA_opt.R")

A <- GSOA_opt(c(1,1,2),s = 3)
cor <- function(a,b){
  return(sum((a-mean(a))*(b-mean(b)))/(sum((a-mean(a))^2)*sum((b-mean(b))^2))^0.5)
}

for (i in 1:7) {
  for(j in (i+1):8)
    print(cor(A[,i],A[,j]))
}
