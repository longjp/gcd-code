## compare continuous causal dantiz with discretized version
rm(list=ls())
set.seed(1234)
library(tidyr)
library(ggplot2)
library(latex2exp)

n <- 100
a0 <- 1
av <- 10
beta <- 1
N <- 1000

est <- matrix(NA_real_,nrow=N,ncol=3)
colnames(est) <- c("GCD","OLS","CD")

GenData <- function(n){
  E <- runif(n)
  h <- rnorm(n)
  x <- 3*h + (av*E + a0)*rnorm(n)
  y <- 9*h + beta*x + rnorm(n)
  return(list(E=E,h=h,x=x,y=y))  
}

for(ii in 1:N){
  dat <- GenData(n)
  E <- dat$E
  x <- dat$x
  y <- dat$y
  Eb <- 1*(E > median(E))
  E <- E - mean(E)
  est[ii,1] <- mean(E*y*x)/mean(E*x^2)
  est[ii,2] <- lm(y~x-1)$coefficients
  num <- mean(x[Eb==1]*y[Eb==1]) - mean(x[Eb==0]*y[Eb==0])
  den <- mean(x[Eb==1]^2) - mean(x[Eb==0]^2)
  est[ii,3] <- num/den
}



n1 <- 100000
dat <- GenData(n1)
E <- dat$E - mean(dat$E)
x <- dat$x
y <- dat$y
## what is asymptotic sd
s <- sqrt((mean(E^2*x^2*(y - x)^2) / mean(E*x^2)^2)/n)


res <- pivot_longer(as.data.frame(est),1:3)
p <- ggplot(res,aes(x=value,fill=name)) +
  geom_histogram(alpha=0.5,position="identity",binwidth=0.1) +
  labs(fill='Estimator') +
  theme(legend.position = c(0.2, 0.7)) +
    labs(x=TeX("$\\widehat{\\beta}$"),
       y=TeX("Count")) +
  coord_cartesian(xlim=c(0,2.5))

pdf("cont_cd.pdf",width=5,height=4)
print(p)
dev.off()

