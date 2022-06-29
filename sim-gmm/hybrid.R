## compare iv, gcd, and hybrid estimator
rm(list=ls())
set.seed(1234)
library(AER)
library(ggplot2)
library(tidyr)
library(latex2exp)
source('funcs.R')

n <- 100
a0 <- 1
beta <- 1
N <- 1000

cnames <- c("GCD","IV","Hybrid")
est <- matrix(NA_real_,nrow=N,ncol=length(cnames))
colnames(est) <- cnames

GenData <- function(n,av,R){
  E <- runif(n)
  h <- rnorm(n)
  x <- h + R*E + (av*E + a0)*rnorm(n)
  y <- h + beta*x + rnorm(n)
  return(list(E=E,h=h,x=x,y=y))  
}

## GCD strong
av <- 5
R <- 1
for(ii in 1:N){
  dat <- GenData(n,av,R)
  E <- dat$E
  x <- dat$x
  y <- dat$y
  E <- as.matrix(E - mean(E),ncol=1)
  x <- as.matrix(x - mean(x),ncol=1)
  y <- y - mean(y)
  est[ii,1] <- ivreg.fit(x,y,E*x)$coefficients
  est[ii,2] <- ivreg.fit(x,y,E)$coefficients
  est[ii,3] <- GCD(E,x,y,center=TRUE,hybrid=TRUE)$coefficients
}


res <- pivot_longer(as.data.frame(est),1:ncol(est))
p1 <- ggplot(res,aes(x=value,fill=name)) +
  geom_histogram(alpha=0.5,position="identity",binwidth=0.1) +
  labs(fill='Estimator') +
  theme(legend.position = c(0.2, 0.7)) +
    labs(x=TeX("$\\widehat{\\beta}$"),
       y=TeX("Count")) +
  coord_cartesian(xlim=c(0,2.5))


pdf("gcd-strong.pdf",width=5,height=4)
print(p1)
dev.off()



## IV strong
av <- 1
R <- 5
for(ii in 1:N){
  dat <- GenData(n,av,R)
  E <- dat$E
  x <- dat$x
  y <- dat$y
  E <- as.matrix(E - mean(E),ncol=1)
  x <- as.matrix(x - mean(x),ncol=1)
  y <- y - mean(y)
  est[ii,1] <- ivreg.fit(x,y,E*x)$coefficients
  est[ii,2] <- ivreg.fit(x,y,E)$coefficients
  est[ii,3] <- GCD(E,x,y,center=TRUE,hybrid=TRUE)$coefficients
}

res <- pivot_longer(as.data.frame(est),1:ncol(est))
p2 <- ggplot(res,aes(x=value,fill=name)) +
  geom_histogram(alpha=0.5,position="identity",binwidth=0.1) +
  labs(fill='Estimator') +
  theme(legend.position = c(0.2, 0.7)) +
    labs(x=TeX("$\\widehat{\\beta}$"),
       y=TeX("Count")) +
  coord_cartesian(xlim=c(0,2.5))

pdf("iv-strong.pdf",width=5,height=4)
print(p2)
dev.off()
