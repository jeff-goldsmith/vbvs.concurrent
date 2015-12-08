#' vb_concurrent
#' 
#' Implements variational bayes for the linear functional concurrent model
#' 
#' @param formula formula for desired regression. should have form \code{ y ~ x1 + x2 + ... + x_k | t}, where \code{t} is the variable that parameterized observed functions
#' @param id.var variable giving subject ID vector
#' @param data optional data frame
#' @param Kt number of spline basis functions for coefficients and FPCs
#' @param Kp number of FPCs to estimate
#' @param v1 tuning parameter; normal slab variance
#' @param standardized logical; are covariates already standardized?
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). if `NULL`, taken to be minium observed value.
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). if `NULL`, taken to be minium observed value.
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' 
#' @importFrom splines bs
#' @importFrom FNN knn.index
#' 
#' @export
#' 
vb_concurrent = function(formula, id.var = NULL, data=NULL, Kt = 5, Kp = 2, v1 = 100,
                         standardized = FALSE, t.min = NULL, t.max = NULL){
  
  ## parse formula
  LHS = as.character(formula)[2]
  RHS = as.character(formula)[3]
  
  if(is.null(id.var)){ stop("Please specify the subject ID variable")}
  if(grep("|", RHS) == 0){ stop("Formula incorrectly specified; needs a '|' to indicate time parameter")}
  
  time.var = strsplit(RHS, "|", fixed = TRUE)[[1]][2] %>% gsub(" ", "", x = ., fixed = TRUE)
  RHS = strsplit(RHS, "|", fixed = TRUE)[[1]][1]
  RHS.mod = paste(RHS, "+", time.var, "+", id.var)
  
  formula.temp = as.formula(paste0(LHS, "~", RHS.mod))
  
  tf <- terms.formula(formula.temp, specials = NULL)
  trmstrings <- attr(tf, "term.labels")
  data.complete = model.frame(tf, data = data)
  
  ## construct theta matrix
  time = data.complete[time.var][,1]
  if(is.null(t.min)) {t.min = min(time)}
  if(is.null(t.max)) {t.max = max(time)}
  
  Theta = t(bs(c(t.min, t.max, time), knots = quantile(time, probs = seq(0, 1, length = Kt - 2))[-c(1,Kt - 2)], 
               intercept=TRUE, degree=3))[,-(1:2)]
  
  formula.model = as.formula(paste0(LHS, "~", RHS))
#  formula.model = as.formula(paste0(LHS, "~ 0+", RHS))
  tf <- terms.formula(formula.model, specials = NULL)
  trmstrings <- attr(tf, "term.labels")
  p = length(trmstrings) + 1
#  p = length(trmstrings)
  mf_fixed = data.model <- model.frame(tf, data = data)
  
  ## normalize variables
  if(!standardized){
    cat("Standardizing Variables \n")
    
    knn.index = as.data.frame(t(knnx.index(time, time, k = round(length(time)*.2))))
    
    for(trm in c(trmstrings)){
#    for(trm in c(LHS, trmstrings)){
      covar.cur = mf_fixed[trm][,1]
      
      mean.fit = sapply(knn.index, function(u){mean(covar.cur[u])})
      sq.resid = (covar.cur - mean.fit)^2
      
      var.fit = sapply(knn.index, function(u){mean(sq.resid[u])})
      
      data.model[trm] = mf_fixed[trm] = (covar.cur - mean.fit) / sqrt(var.fit)
      cat(".")
    }
    cat("\n")
  }
  
  data.model[id.var] = data.complete[id.var]
  data.model[time.var] = data.complete[time.var]
  
  Y = data.model[LHS][,1]
  subj.id = data.model[id.var][,1]
  subjs = unique(subj.id)
  I = length(subjs)
  J = dim(data.model)[1]
  
  ## construct Xstar
  cat("Constructing Xstar; doing data organization \n")
  
  Xstar = tXsXs = vector("list", length = I)
  XStXS = matrix(0, Kt*p, Kt*p)
  XStY = matrix(0, Kt*p, 1)
  sumtXsXs = matrix(0, Kt*p, Kt*p)
  one.p = matrix(1, p, 1)
  one.Kt = matrix(1, Kt, 1)
  for(i in 1:I){
    index = which(subj.id == subjs[i])
    
    X.cur = t(model.matrix(mf_fixed, data = slice(data.model, index)))
    Theta.cur = Theta[,index]
    Y.cur = Y[index]
    XS.cur = Xstar[[i]] = t(kronecker(X.cur, one.Kt) * kronecker(one.p, Theta.cur))
    XStXS = XStXS + crossprod(XS.cur)
    XStY = XStY + t(XS.cur) %*% as.matrix(Y.cur)
    sumtXsXs = sumtXsXs + t(XS.cur)%*% XS.cur
  }
  
  ## hyper parameters for inverse gaussians, bernoulli
  A = .5
  B = .5

  ## matrices to to approximate paramater values
  sigma.q.BW = vector("list", p)
  for(k in 1:p){
    sigma.q.BW[[k]] = diag(1, Kt)
  }
  mu.q.BW = matrix(0, nrow = Kt, ncol = p)  
  
  mu.q.gamma = rep(1, p)  
  mu.q.dinv = kronecker(diag(mu.q.gamma/v1, p, p), diag(1, Kt, Kt))
  
  sigma.q.Bpsi = vector("list", Kp)
  for(k in 1:Kp){
    sigma.q.Bpsi[[k]] = diag(1, Kt)
  }
  mu.q.Bpsi = matrix(0, nrow = Kt, ncol = Kp)  
  
  sigma.q.C = vector("list", I)
  for(k in 1:I){
    sigma.q.C[[k]] = diag(1, Kp)
  }
  set.seed(1)
  mu.q.C = matrix(rnorm(I*Kp, 0, 1), I, Kp)
  
  b.q.lambda.Bpsi = rep(1, Kp)
  b.q.sigma.me = 1

  ## initialize estimates of fixed and pca effects
  fixef.est = rep(0, J)
  pcaef.est = rep(0, J)
  
  cat("Beginning Algorithm \n")

  for(iter in 1:10){
    
    ###############################################################
    ## update regression coefficients
    ###############################################################
  
    XStY = matrix(0, Kt*p, 1)
    for(i in 1:I){
      index = which(subj.id == subjs[i])
      Y.cur = Y[index]
      pcaef.cur = pcaef.est[index]
      XS.cur = Xstar[[i]]
      XStY = XStY + t(XS.cur) %*% as.matrix(Y.cur - pcaef.cur)
    }
    
    sigma.q.beta = solve(as.numeric((A + J/2)/(b.q.sigma.me)) * XStXS + mu.q.dinv)
    mu.q.beta = matrix(sigma.q.beta %*% (as.numeric((A + J/2)/(b.q.sigma.me)) * XStY), nrow = Kt, ncol = p)
    
    beta.cur = t(mu.q.beta) %*% (Theta)
    
    for(i in 1:I){
      index = which(subj.id == subjs[i])
      fixef.est[index] = Xstar[[i]] %*% as.vector(mu.q.beta)
      # equivalent formulation: apply(X[[subj]] * beta.cur, 2, sum)
    }
    
    ###############################################################
    ## update gammas -- not in this version!
    ###############################################################
  

    ###############################################################
    ## update b-spline parameters for PC basis functions
    ###############################################################
    
    sigma.q.Bpsi = solve( 
      kronecker(diag(1, Kt, Kt), diag((A+Kt/2)/b.q.lambda.Bpsi)) + 
        as.numeric((A + J/2)/(b.q.sigma.me)) * f_sum(mu.q.c = mu.q.C, sig.q.c = sigma.q.C, theta = Theta, subj.id = subj.id)
    )
    mu.q.Bpsi = matrix(((A + J/2)/(b.q.sigma.me)) * sigma.q.Bpsi %*% f_sum2(y = Y, fixef = fixef.est, subj.id = subj.id, mu.q.c = mu.q.C, kt = Kt, theta = Theta), nrow = Kt, ncol = Kp)
    
    psi.cur = t(mu.q.Bpsi) %*% (Theta)

    ###############################################################
    ## scores for each individual
    ###############################################################
    
    for(i in 1:I){
      index = which(subj.id == subjs[i])
      
      Theta_i = Theta[,index]
      sigma.q.C[[i]] = solve( 
        diag(1, Kp, Kp ) +
          ((A + J/2)/(b.q.sigma.me)) * (f_trace(Theta_i = Theta_i, Sig_q_Bpsi = sigma.q.Bpsi, Kp = Kp, Kt = Kt) + 
                                          t(mu.q.Bpsi) %*% Theta_i %*% t(Theta_i) %*% mu.q.Bpsi)
      )
      
      mu.q.C[i,] = ((A + J/2)/(b.q.sigma.me)) * sigma.q.C[[i]] %*% as.matrix(psi.cur[,index]) %*%  (Y[index] - fixef.est[index] )
    }
    
    for(i in 1:I){
      index = which(subj.id == subjs[i])
      pcaef.est[index] = mu.q.C[i,] %*% psi.cur[,index]
    }
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    ## measurement error variance
    resid = Y - fixef.est - pcaef.est
    b.q.sigma.me = as.numeric(B + .5 * (crossprod(resid) + 
                                          sum(diag(sumtXsXs %*% sigma.q.beta)) + 
                                          f_sum4(mu.q.c= mu.q.C, sig.q.c = sigma.q.C, mu.q.bpsi = mu.q.Bpsi, sig.q.bpsi = sigma.q.Bpsi, theta = Theta, subj.id = subj.id)) )
    
    
    ## lambda for FPCA basis functions
    for(K in 1:Kp){
      b.q.lambda.Bpsi[K] = B + .5 * (t(mu.q.Bpsi[,K]) %*% diag(1, Kt, Kt) %*% mu.q.Bpsi[,K] + 
                                       sum(diag(diag(1, Kt, Kt) %*% sigma.q.Bpsi[(Kt*(K-1)+1):(Kt*K),(Kt*(K-1)+1):(Kt*K)])))
    }
    
    cat(".")
    
  }
  
  cat("\n")
  

  ## get coefficient functions over a common grid
  rownames(beta.cur) = c("int", trmstrings)
#  rownames(beta.cur) = c(trmstrings)
  beta.cur = t(beta.cur) %>% as.data.frame() %>%
    mutate(t = time) %>% 
    arrange(t) %>% unique() %>%
    subset(select= c("t", "int", trmstrings))
#    subset(select= c("t", trmstrings))
  
  ## export fitted values
  Yhat = fixef.est + pcaef.est 
  
  ## export variance components
  sigeps.pm = 1 / as.numeric((A + J/2)/(b.q.sigma.me))
  
  ## do svd to get rotated fpca basis
  rownames(psi.cur) = paste0("Psi.", 1:Kp)
  psi.cur = t(psi.cur) %>% as.data.frame() %>%
    mutate(t = time) %>% 
    arrange(t) %>% unique() %>%
    subset(select= c("t", paste0("Psi.", 1:Kp)))
  
  temp = svd(psi.cur[,-1])
  psi.cur[,-1] = temp$u
  lambda.pm = temp$d
  
  ret = list(beta.cur, psi.cur, mu.q.beta, lambda.pm, Yhat, sigeps.pm, pcaef.est, fixef.est, 
             data.complete, data.model, formula, formula.model, time.var, id.var, Kt, Kp, standardized,
             t.min, t.max)
  names(ret) = c("beta.pm", "psi.pm", "spline.coef.est", "lambda.pm", "Yhat", "sigeps.pm", "pcaef", "fixef",
                 "data", "data.model", "formula", "formula.model", "time.var", "id.var", "Kt", "Kp", "standardized",
                 "t.min", "t.max")
  
  class(ret) = "concurFLM"
  
  ret
  
}

