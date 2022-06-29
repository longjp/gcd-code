## functions for fitting the GCD and hybrid estimator

## fits gcd given kronecker product and weight matrix
GCD_fit <- function(XkrE,X,y,What){
  A <- t(X)%*%XkrE%*%What%*%t(XkrE)%*%X
  B <- colSums(t(t(X)%*%XkrE%*%What%*%t(XkrE))*y)
  solve(A,B)
}

## computes rowwise kronecker product of A and B
RowwiseKron <- function(A,B){
  return(A[,rep(seq(ncol(A)),each = ncol(B)),drop=FALSE]*B[,rep(seq(ncol(B)),ncol(A)),drop=FALSE])
}

## fits GCD or hybrid (CD + IV constraints) estimator
##
## note: in overidentified case (where number of constraints
## exceeds number of predictors), the estimator implements a
## two step algorithm for estimating the optimal weight matrix.
## this two step estimator uses a two stage least squares like
## initial estimate of the weight matrix W
##
## returns point estimates and covariance matrix
##
##
## arguments
##       E : n x q matrix containing instruments/environments
##       X : n x p matrix of predictors
##       y : length n vector with response
##  center : should y and X columns be set to mean 0
##  hybrid : should instrumental variable constraints be added to estimator
##
##  note on coding E: if using r environments 0, 1, . . ., r-1
##                    column j of E should be 1 if
##                    observation is environment j, 0 o.w.
GCD <- function(E,X,y,center=TRUE,hybrid=FALSE){
  ## normalize data
  n <- nrow(E)
  E <- t(t(E) - colMeans(E))
  if(center){
    X <- t(t(X) - colMeans(X))
    y <- y - mean(y)
  }
  ## compute X dot E (kronecker product of X and E)
  XkrE <- RowwiseKron(X,E)
  ## for hybrid estimator, IV orthogonality constraints are added
  if(hybrid){
    XkrE <- cbind(XkrE,E)
  }
  ## initial weight matrix, motivated by tsls
  What <- solve((1/n)*t(XkrE)%*%XkrE)
  ## in overidentified case, need two steps
  if(ncol(XkrE)>ncol(X)){
    betahat <- GCD_fit(XkrE,X,y,What)
    deltahat <- y - colSums(t(X)*betahat)
    temp <- XkrE*deltahat
    What <- solve((1/n)*t(temp)%*%temp)
  }
  ## compute estimate and estimate of asymptotic variance
  betahat <- GCD_fit(XkrE,X,y,What)
  deltahat <- y - colSums(t(X)*betahat)
  temp <- XkrE*deltahat
  V <- t(temp)%*%temp/n
  M <- -(1/n)*t(XkrE)%*%X
  betahat_var <- solve(t(M)%*%solve(V)%*%M)/n
  return(list(coefficients=betahat,
              cov=betahat_var))
}
