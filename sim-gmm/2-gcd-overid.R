## fit GCD
rm(list=ls())
set.seed(1234)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(kableExtra)
source("funcs.R")

n <- 200

sig <- 3
del <- 5
N <- 500


methods <- c("GCD","GCDE1","GCDE2","OLS")
res <- vector("list",length=length(methods))
names(res) <- methods
for(ii in 1:length(res)){
  res[[ii]] <- vector("list",N)
}
for(ii in 1:N){
  E1 <- rbinom(n,size=1,prob=1/2)
  E2 <- runif(n)
  eta <- matrix(rnorm(n*5),ncol=5)
  x2 <- eta[,1] + (1+sig*E1+del*E2)*eta[,4]
  y <- x2 + eta[,1] + eta[,2]
  x1 <- y + x2 + (1+sig*E1)*eta[,3]
  x3 <- x1 + eta[,1] + (1+del*E2)*eta[,5]
  E <- cbind(E1,E2)
  X <- cbind(x1,x2,x3)
  res$GCD[[ii]] <- GCD(E,X,y)
  res$GCDE1[[ii]] <- GCD(as.matrix(E1,ncol=1),X,y)
  res$GCDE2[[ii]] <- GCD(as.matrix(E2,ncol=1),X,y)
  a <- lm(y~X-1)
  res$OLS[[ii]] <- list(coefficients=a$coefficients,cov=vcov(a)) 
}



out <- matrix(NA_real_,nrow=4,ncol=3)
rownames(out) <- names(res)
for(ii in 1:length(res)){
  temp <- matrix(NA,ncol=3,nrow=length(res[[ii]]))
  for(jj in 1:length(res[[ii]])){
    f <- res[[ii]][[jj]]
    coef <- f$coefficients
    temp[jj,] <- abs((coef-c(0,1,0))/sqrt(diag(f$cov)))<2
  }
  out[ii,] <- colMeans(temp)
}
out <- round(out,2)
colnames(out) <- c("$\\beta_1$","$\\beta_2$","$\\beta_3$")



wid <- matrix(NA_real_,nrow=4,ncol=3)
rownames(wid) <- names(res)
for(ii in 1:length(res)){
  temp <- matrix(NA,ncol=3,nrow=length(res[[ii]]))
  for(jj in 1:length(res[[ii]])){
    f <- res[[ii]][[jj]]
    temp[jj,] <- 4*sqrt(diag(f$cov))
  }
  wid[ii,] <- apply(temp,2,median)
}
wid <- round(wid,2)
colnames(wid) <- c("$\\beta_1$","$\\beta_2$","$\\beta_3$")



out <- as.data.frame(cbind(out,wid))

a <- kbl(out,caption="Coverage and Median Width of Confidence Intervals for Different Estimators \\label{tab:overid}",
         format="latex",escape=FALSE) %>%
  add_header_above(c(" "=1,"Coverage"=3,"Median Width"=3))
save_kable(a,file="gcd-overid.tex")

# 
# 
# a <- solve(f$cov)
# a.eig <- eigen(a)
# a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
# (a.sqrt%*%matrix(coef0,ncol=1))[,1]
# 
# 
# 
# fit <- GCD(E,X,y)
# colnames(res) <- paste0("beta",1:3)
# res <- as.data.frame(res)
# 
# 
# ## x1 versus x2
# elCD <- ellipse::ellipse(fit$cov[1:2,1:2],centre=c(0,1))
# elCD <- as.data.frame(elCD)
# colnames(elCD) <- c("x1","x2")
# xlim <- range(elCD$x1)
# ylim <- range(elCD$x2)
# p1 <- ggplot(res,aes(x=beta1,y=beta2)) +
#   geom_point(shape=21,alpha=0.5) +
#   geom_hline(yintercept=1,alpha=0.3) + 
#   geom_vline(xintercept=0,alpha=0.3) + 
#   coord_cartesian(xlim=xlim,ylim=ylim) +
#   geom_path(data=elCD,mapping=aes(x=x1,y=x2),inherit.aes=FALSE) +
#   labs(x=TeX("$\\widehat{\\beta}_1$"),y=TeX("$\\widehat{\\beta}_2$"))
# print(p1)
# 
# 
# 
# 
# ## x1 versus x3
# elCD <- ellipse::ellipse(fit$cov[c(1,3),c(1,3)],centre=c(0,0))
# elCD <- as.data.frame(elCD)
# colnames(elCD) <- c("x1","x3")
# xlim <- range(elCD$x1)
# ylim <- range(elCD$x3)
# p2 <- ggplot(res,aes(x=beta1,y=beta3)) +
#   geom_point(shape=21,alpha=0.5) +
#   geom_hline(yintercept=0,alpha=0.3) + 
#   geom_vline(xintercept=0,alpha=0.3) + 
#   coord_cartesian(xlim=xlim,ylim=ylim) +
#   geom_path(data=elCD,mapping=aes(x=x1,y=x3),inherit.aes=FALSE) +
#   labs(x=TeX("$\\widehat{\\beta}_1$"),y=TeX("$\\widehat{\\beta}_3$"))
# print(p2)
# 
# 
# 
# pdf("gcd-overid-1.pdf",width=4,height=4)
# plot(p1)
# dev.off()
# 
# pdf("gcd-overid-2.pdf",width=4,height=4)
# plot(p2)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## analyze GCD results
# res <- pivot_longer(as.data.frame(res),1:ncol(res))
# p <- ggplot(res,aes(x=value,fill=name)) +
#   geom_histogram(alpha=0.5,position="identity",binwidth=0.1) +
#   labs(fill='Estimator') +
#   theme(legend.position = c(0.2, 0.7)) +
#     labs(x=TeX("$\\widehat{\\beta}$"),
#        y=TeX("Count")) +
#   coord_cartesian(xlim=c(-1,2))
# 
# print(p)
# 
# 
# 
# ## analyze OLS results
# res_lm <- pivot_longer(as.data.frame(res_lm),1:ncol(res_lm))
# p <- ggplot(res_lm,aes(x=value,fill=name)) +
#   geom_histogram(alpha=0.5,position="identity",binwidth=0.1) +
#   labs(fill='Estimator') +
#   theme(legend.position = c(0.2, 0.7)) +
#     labs(x=TeX("$\\widehat{\\beta}$"),
#        y=TeX("Count")) +
#   coord_cartesian(xlim=c(-1,2))
# 
# print(p)
# 
# 
# 
# ## analyze OLS results
# res_E1 <- pivot_longer(as.data.frame(res_E1),1:ncol(res_E1))
# p <- ggplot(res_E1,aes(x=value,fill=name)) +
#   geom_histogram(alpha=0.5,position="identity",binwidth=0.1) +
#   labs(fill='Estimator') +
#   theme(legend.position = c(0.2, 0.7)) +
#     labs(x=TeX("$\\widehat{\\beta}$"),
#        y=TeX("Count")) +
#   coord_cartesian(xlim=c(-1,2))
# 
# print(p)
# 
