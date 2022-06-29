## tests of funcs.R
## TODO
## 1) compare GCD with CD in InvariantCausalPrediction Package
##    on simulation in paper
## 2) compare with own code
## 3) for continuous and hybrid simulations, use funcs.R 
##    implementation of GCD. make sure we get same results as now
## 4) relationship between cDantzig function below and code in ICP package, are they different
## 5) change section GenData function so that GCD is more stable
## 6) why do GCD and cDantzig have different asymptotic variances for simulation
rm(list=ls())
source('funcs.R')
library(ggplot2)
set.seed(1234)

#### older implementation of CD
cDantzig <- function(X,Y,ExpInd){
  Y <- Y - mean(Y)
  X <- t(t(X) - colMeans(X))
  X1 <- X[ExpInd == 1,,drop=FALSE]
  Y1 <- Y[ExpInd == 1]
  X0 <- X[ExpInd == 0,,drop=FALSE]
  Y0 <- Y[ExpInd == 0]
  Ginv <- solve(t(X1) %*% X1/nrow(X1) - t(X0) %*% X0/nrow(X0))
  Z <- t(X1) %*% Y1/nrow(X1) - t(X0) %*% Y0/nrow(X0)
  beta <- Ginv%*%Z
  v00 <- cov(X0*as.vector(Y0-X0%*%beta)) / nrow(X0)
  v11 <- cov(X1*as.vector(Y1-X1%*%beta)) / nrow(X1)
  Va <- Ginv%*%(v00 + v11)%*%Ginv
  # v0 <- X0*as.vector(Y0-X0%*%beta)
  # v0 <- t(v0)%*%v0/nrow(X0)^2
  # v1 <- X1*as.vector(Y1-X1%*%beta)
  # v1 <- t(v1)%*%v1/nrow(X1)^2
  # Va <- Ginv%*%(v0+v1)%*%Ginv
  return(list(beta=beta,V=Va))
}


## COMP 1
## on 1-d problem compare
## GCD function with simple 1-d formula
GenData <- function(n){
  e <- runif(n)
  h <- rnorm(n)
  x <- 3*h + (av*e + a0)*rnorm(n)
  y <- 9*h + beta*x + rnorm(n)
  return(list(e=e,h=h,x=x,y=y))  
}

n <- 1000
a0 <- 1
av <- 10
beta <- 1


dat <- GenData(n)
e <- dat$e
x <- dat$x
y <- dat$y

X <- matrix(x,ncol=1)
E <- matrix(e,ncol=1)
GCD(E,X,y)
e <- e - mean(e)
x <- x - mean(x)
y <- y - mean(y)
mean(e*y*x)/mean(e*x^2)


N <- 5000
res <- rep(NA_real_,N)
for(ii in 1:N){
  dat <- GenData(n)
  e <- dat$e
  x <- dat$x
  y <- dat$y
  X <- matrix(x,ncol=1)
  E <- matrix(e,ncol=1)
  fit <- GCD(E,X,y)
  res[ii] <- (fit$coefficients-beta)/sqrt(fit$cov)
}
df <- data.frame(res)
# histogram with normal density
ggplot(df, aes(x = res)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = dnorm)



## COMP 2
GenData <- function(n){
  e1 <- rnorm(n)
  e2 <- rnorm(n)
  ExpInd <- rbinom(n,size=1,prob=1/2)
  e1 <- e1*ExpInd
  e2 <- e2*ExpInd
  h <- rnorm(n)
  x1 <- e1*rnorm(n,sd=3) + e2*rnorm(n,sd=6) + rnorm(n) + h
  x2 <- e1*rnorm(n,sd=2) + e2*rnorm(n,sd=8) + rnorm(n) + 2*h
  y <- x1 + 2*h + rnorm(n)
  return(list(y=y,x1=x1,x2=x2,e1=e1,e2=e2,ExpInd=ExpInd))
}

n <- 1000 ## number of observations
dat <- GenData(n)
X <- cbind(dat$x1,dat$x2)
y <- dat$y
f1 <- cDantzig(X,y,dat$ExpInd)
E <- matrix(dat$ExpInd,ncol=1)
f2 <- GCD(E,X,y)

## should match
f1$beta
f2$coefficients





n <- 1000 ## number of observations
N <- 5000
res <- matrix(NA_real_,nrow=N,ncol=2)
for(ii in 1:N){
  dat <- GenData(n)
  X <- cbind(dat$x1,dat$x2)
  y <- dat$y
  E <- matrix(dat$ExpInd,ncol=1)
  f2 <- GCD(E,X,y)
  coef0 <- f2$coefficients - c(1,0)
  a <- solve(f2$cov)
  a.eig <- eigen(a)
  a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
  res[ii,] <- (a.sqrt%*%matrix(coef0,ncol=1))[,1]
}

## should be 2d standard normal, looks good
plot(res[,1],res[,2])

## check that normalized sum follows standard normal
res <- (res[,1] + res[,2])/sqrt(2)
df <- data.frame(res)
# Histogram with kernel density
ggplot(df, aes(x = res)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = dnorm)

