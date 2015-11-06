#' ols_concurrent
#' 
#' Implements ordinary least squares for the linear functional concurrent model
#' 
#' @param Y matrix of functional responses
#' @param X list of matrices, each of which contains the functional predictor for an individual subject
#' @param Kt number of spline basis functions for coefficients and FPCs
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' @export
ols_concurrent = function(Y, X, Kt=5){
  
  ## subject covariates
  I = dim(Y)[1]
  D = dim(Y)[2]
  p = dim(X[[1]])[1]
  
  ## bspline basis and penalty matrix
  Theta = t(bs(1:D, df=Kt, intercept=TRUE, degree=3))
  
  ## define Xstar
  Xstar = vector("list", length = I)
  one.p = matrix(1, p, 1)
  one.Kt = matrix(1, Kt, 1)
  for(i in 1:I){
    Xstar[[i]] = t(kronecker(X[[i]], one.Kt) * kronecker(one.p, Theta))
  }
  
  ## data organization; these computations only need to be done once
  XStXS = matrix(0, Kt*p, Kt*p)
  XStY = matrix(0, Kt*p, 1)
  sumtXsXs = matrix(0, Kt*p, Kt*p)
  for(i in 1:I){
    obs.points = which(!is.na(Y[i, ]))
    XS.cur = Xstar[[i]][obs.points,]
    XStXS = XStXS + crossprod(XS.cur)
    XStY = XStY + t(XS.cur) %*% as.matrix(Y[i,obs.points])
  }
  
  bhat = matrix(solve(XStXS)%*%XStY, Kt, p)
  betahat = t(bhat) %*% (Theta)
  
  ret = list(betahat)
  names(ret) = c("betahat")
  
  ret
  
}