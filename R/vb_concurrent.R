#' vb_concurrent
#' 
#' Implements variational bayes for the linear functional concurrent model
#' 
#' @param Y matrix of functional responses
#' @param X list of matrices, each of which contains the functional predictor for an individual subject
#' @param Kt number of spline basis functions for coefficients and FPCs
#' @param Kp number of FPCs to estimate
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' @export
vb_concurrent = function(Y, X, Kt = 5, Kp = 2){
  
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
  
  ## hyper parameters for inverse gaussians, bernoulli
  v0 = .01
  v1 = 100
  A = .5
  B = .5
  Atheta = .5
  Btheta = .5
  
  
  ## matrices to to approximate paramater values
  sigma.q.BW = vector("list", p)
  for(k in 1:p){
    sigma.q.BW[[k]] = diag(1, Kt)
  }
  mu.q.BW = matrix(0, nrow = Kt, ncol = p)  
  
  mu.q.gamma = rep(1, p)  
  mu.q.dinv = kronecker(diag((1-mu.q.gamma)/v0 + mu.q.gamma/v1, p, p), diag(1, Kt, Kt))
  
  sigma.q.Bpsi = vector("list", Kp)
  for(k in 1:Kp){
    sigma.q.Bpsi[[k]] = diag(1, Kt)
  }
  mu.q.Bpsi = matrix(0, nrow = Kt, ncol = Kp)  
  
  sigma.q.C = vector("list", I)
  for(k in 1:I){
    sigma.q.C[[k]] = diag(1, Kp)
  }
  mu.q.C = matrix(rnorm(I*Kp, 0, 1), I, Kp)
  
  mu.q.ltheta = digamma(Atheta + sum(mu.q.gamma)) - digamma(Atheta + Btheta + p)
  mu.q.l1.theta = digamma(Btheta + p - sum(mu.q.gamma)) - digamma(Atheta + Btheta + p)
  
  b.q.lambda.Bpsi = rep(1, Kp)
  b.q.sigma.me = 1
  
  ## data organization; these computations only need to be done once
  XStXS = matrix(0, Kt*p, Kt*p)
  XStY = matrix(0, Kt*p, 1)
  sumtXsXs = matrix(0, Kt*p, Kt*p)
  for(i in 1:I){
    obs.points = which(!is.na(Y[i, ]))
    XS.cur = Xstar[[i]][obs.points,]
    XStXS = XStXS + crossprod(XS.cur)
    XStY = XStY + t(XS.cur) %*% as.matrix(Y[i,obs.points])
    sumtXsXs = sumtXsXs + t(XS.cur)%*% XS.cur
  }
  
  
  Y.vec = as.vector(t(Y))
  obspts.vec = !is.na(Y.vec)
  Y.vec = Y.vec[obspts.vec]
  J = sum(obspts.vec)
  
  
  ## initialize estimates of fixed and pca effects
  fixef.cur = matrix(0, nrow = I, ncol = D)
  pcaef.cur = matrix(0, I, D)
  
  for(iter in 1:10){
    
    ###############################################################
    ## update regression coefficients
    ###############################################################
  
    XStY = matrix(0, Kt*p, 1)
    for(subj in 1:I){
      obs.points = which(!is.na(Y[subj, ]))
      XS.cur = Xstar[[subj]][obs.points,]
      XStY = XStY + t(XS.cur) %*% as.matrix(Y[subj,obs.points] - pcaef.cur[subj,obs.points])
    }
    
    sigma.q.beta = solve(as.numeric((A + J/2)/(b.q.sigma.me)) * XStXS + mu.q.dinv)
    mu.q.beta = matrix(sigma.q.beta %*% (as.numeric((A + J/2)/(b.q.sigma.me)) * XStY), nrow = Kt, ncol = p)
    
    beta.cur = t(mu.q.beta) %*% (Theta)
    
    for(subj in 1:I){
      fixef.cur[subj,] = Xstar[[subj]] %*% as.vector(mu.q.beta)
          # equivalent formulation: apply(X[[subj]] * beta.cur, 2, sum)
    }
    
    ###############################################################
    ## update gammas -- not in this version!
    ###############################################################
  

    ###############################################################
    ## update b-spline parameters for PC basis functions
    ###############################################################
    
    mean.cur = as.vector(t(fixef.cur))[obspts.vec]
    designmat = kronecker(mu.q.C, t(Theta))[obspts.vec,]
    
    sigma.q.Bpsi = solve( 
        kronecker(diag(1, Kt, Kt), diag((A+Kt/2)/b.q.lambda.Bpsi)) + 
        as.numeric((A + J/2)/(b.q.sigma.me)) * f.sum(mu.q.c = mu.q.C, sig.q.c = sigma.q.C, theta = Theta, obspts.mat = !is.na(Y))
      )
    mu.q.Bpsi = matrix(((A + J/2)/(b.q.sigma.me)) * sigma.q.Bpsi %*% f.sum2(y = Y, fixef = fixef.cur, mu.q.c = mu.q.C, kt = Kt, theta = Theta), nrow = Kt, ncol = Kp)
    
    psi.cur = t(mu.q.Bpsi) %*% (Theta)
    ppT = (psi.cur) %*% t(psi.cur)

    ###############################################################
    ## scores for each individual
    ###############################################################
    
    for(subj in 1:I){
      obs.points = which(!is.na(Y[subj, ]))
      Theta_i = Theta[,obs.points]
      sigma.q.C[[subj]] = solve( 
        diag(1, Kp, Kp ) +
        ((A + J/2)/(b.q.sigma.me)) * (f.trace(Theta_i = Theta_i, Sig_q_Bpsi = sigma.q.Bpsi, Kp = Kp, Kt = Kt) + 
                                      t(mu.q.Bpsi) %*% Theta_i %*% t(Theta_i) %*% mu.q.Bpsi)
      )
      
      mu.q.C[subj,] = ((A + J/2)/(b.q.sigma.me)) * sigma.q.C[[subj]] %*% as.matrix(psi.cur[,obs.points]) %*%  (Y[subj,obs.points] - fixef.cur[subj,obs.points] )
    }
    
    pcaef.cur =  as.matrix(mu.q.C %*% psi.cur)
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    ## measurement error variance
    resid = as.vector(Y - fixef.cur - pcaef.cur)
    b.q.sigma.me = as.numeric(B + .5 * (crossprod(resid[!is.na(resid)]) + 
                                          sum(diag(sumtXsXs %*% sigma.q.beta)) + 
                                          f.sum4(mu.q.c= mu.q.C, sig.q.c = sigma.q.C, mu.q.bpsi = mu.q.Bpsi, sig.q.bphi = sigma.q.Bpsi, theta= Theta, obspts.mat = !is.na(Y))) )
    
    
    ## lambda for FPCA basis functions
    for(K in 1:Kp){
      b.q.lambda.Bpsi[K] = B + .5 * (t(mu.q.Bpsi[,K]) %*% diag(1, Kt, Kt) %*% mu.q.Bpsi[,K] + 
                                sum(diag(diag(1, Kt, Kt) %*% sigma.q.Bpsi[(Kt*(K-1)+1):(Kt*K),(Kt*(K-1)+1):(Kt*K)])))
    }
    
#    cat(iter, "\n")
  }

  ## export fitted values
  Yhat = fixef.cur + pcaef.cur

  ## export variance components
  sigeps.pm = 1 / as.numeric((A + J/2)/(b.q.sigma.me))
  
  ## do svd to get rotated fpca basis
  temp = svd(t(psi.cur))
  psi.cur = t(temp$u)
  lambda.pm = temp$d
  
  ret = list(beta.cur, mu.q.gamma, psi.cur, lambda.pm, Yhat, sigeps.pm, pcaef.cur)
  names(ret) = c("beta.pm", "gamma.pm", "psi.pm", "lambda.pm", "Yhat", "sigeps.pm", "pcaef")
  
  ret
  
}