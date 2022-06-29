## CD is consistent as long as skew is non-zero
rm(list=ls())

n <- 100000
E <- sample(c(-1/9,1),size=n,prob=c(9/10,1/10),replace=TRUE)
N <- 1000
res <- matrix(0,ncol=2,nrow=N)
for(ii in 1:N){
  h <- rnorm(n)
  X <- E + h + rnorm(n)
  Y <- X + h + rnorm(n)
  res[ii,1] <- mean(Y*E) / mean(X*E)
  res[ii,2] <- mean(Y*X*E) / mean(X^2*E)
}


summary(res[,1])
summary(res[,2])
