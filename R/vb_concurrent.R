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
#' @importFrom stats model.frame
#' @importFrom splines bs
#' 
#' @export
#' 
vb_concurrent = function(formula, id.var = NULL, data=NULL, Kt = 5, Kp = 2, v1 = 100,
                         standardized = FALSE, t.min = NULL, t.max = NULL){
  
  ## parse formula
  LHS = as.character(formula)[2]
  RHS = as.character(formula)[3]
  
  if (is.null(id.var)) { stop("Please specify the subject ID variable")}
  if (grep("|", RHS) == 0) { stop("Formula incorrectly specified; needs a '|' to indicate time parameter")}
  
  time.var = strsplit(RHS, "|", fixed = TRUE)[[1]][2] %>% gsub(" ", "", x = ., fixed = TRUE)
  RHS = strsplit(RHS, "|", fixed = TRUE)[[1]][1]
  RHS.mod = paste(RHS, "+", time.var, "+", id.var)
  
  formula.temp = as.formula(paste0(LHS, "~", RHS.mod))
  
  tf <- terms.formula(formula.temp, specials = NULL)
  trmstrings <- attr(tf, "term.labels")
  data.complete = model.frame(tf, data = data)
  
  ## construct theta matrix
  time = data.complete[time.var][,1]
  if (is.null(t.min)) {t.min = min(time)}
  if (is.null(t.max)) {t.max = max(time)}
  
  Theta = t(bs(c(t.min, t.max, time), knots = quantile(time, probs = seq(0, 1, length = Kt - 2))[-c(1,Kt - 2)], 
               intercept = TRUE, degree = 3))[,-(1:2)]
  
#  formula.model = as.formula(paste0(LHS, "~", RHS))
  formula.model = as.formula(paste0(LHS, "~ 0+", RHS))
  tf <- terms.formula(formula.model, specials = NULL)
  trmstrings <- attr(tf, "term.labels")
#  p = length(trmstrings) + 1
  p = length(trmstrings)
  mf_fixed = data.model = model.frame(tf, data = data)
  mt_fixed = attr(mf_fixed, "terms")
  
  data.model[id.var] = data.complete[id.var]
  data.model[time.var] = data.complete[time.var]
  
  ## normalize variables
  if (!standardized) {
    
    cat("Standardizing Variables \n")
    data.model = standardize_variables(data.original = data.model, time.var = time.var, trmstrings = trmstrings, LHS = LHS)
    cat("\n")
    
  }
  
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
  for (i in 1:I) {
    index = which(subj.id == subjs[i])
    
    X.cur = t(model.matrix(mf_fixed, data = slice(data.model, index)))
    Theta.cur = Theta[,index]
    Y.cur = Y[index]
    XS.cur = Xstar[[i]] = t(kronecker(X.cur, one.Kt) * kronecker(one.p, Theta.cur))
    XStXS = XStXS + crossprod(XS.cur)
    XStY = XStY + t(XS.cur) %*% as.matrix(Y.cur)
    sumtXsXs = sumtXsXs + t(XS.cur) %*% XS.cur
  }
  
  ## hyper parameters for inverse gaussians, bernoulli
  A = .5
  B = .5

  ## matrices to to approximate paramater values
  mu.q.gamma = rep(1, p)  
  mu.q.dinv = kronecker(diag(mu.q.gamma/v1, p, p), diag(1, Kt, Kt))
  mu.q.Bpsi = matrix(0, nrow = Kt, ncol = Kp)  
  
  sigma.q.C = vector("list", I)
  for (k in 1:I) {
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

  for (iter in 1:10) {
    
    ###############################################################
    ## update regression coefficients
    ###############################################################
  
    XStY = matrix(0, Kt*p, 1)
    for (i in 1:I) {
      index = which(subj.id == subjs[i])
      Y.cur = Y[index]
      pcaef.cur = pcaef.est[index]
      XS.cur = Xstar[[i]]
      XStY = XStY + t(XS.cur) %*% as.matrix(Y.cur - pcaef.cur)
    }
    
    sigma.q.beta = solve(as.numeric((A + J/2)/(b.q.sigma.me)) * XStXS + mu.q.dinv)
    mu.q.beta = matrix(sigma.q.beta %*% (as.numeric((A + J/2)/(b.q.sigma.me)) * XStY), nrow = Kt, ncol = p)
    
    beta.cur = t(mu.q.beta) %*% (Theta)
    
    for (i in 1:I) {
      index = which(subj.id == subjs[i])
      fixef.est[index] = Xstar[[i]] %*% as.vector(mu.q.beta)
      # equivalent formulation: apply(X[[subj]] * beta.cur, 2, sum)
    }
    
    ###############################################################
    ## update gammas -- not in this version!
    ###############################################################
  
    if (Kp > 0) {
      
      ###############################################################
      ## update b-spline parameters for PC basis functions
      ###############################################################
      
      sigma.q.Bpsi = solve( 
        kronecker(diag(1, Kt, Kt), diag((A + Kt/2)/b.q.lambda.Bpsi, Kp, Kp)) + 
          as.numeric((A + J/2)/(b.q.sigma.me)) * f_sum(mu.q.c = mu.q.C, sig.q.c = sigma.q.C, theta = Theta, subj.id = subj.id)
      )
      mu.q.Bpsi = matrix(((A + J/2)/(b.q.sigma.me)) * sigma.q.Bpsi %*% f_sum2(y = Y, fixef = fixef.est, subj.id = subj.id, mu.q.c = mu.q.C, kt = Kt, theta = Theta), nrow = Kt, ncol = Kp)
      
      psi.cur = t(mu.q.Bpsi) %*% (Theta)
      
      ###############################################################
      ## scores for each individual
      ###############################################################
      
      for (i in 1:I) {
        index = which(subj.id == subjs[i])
        
        Theta_i = Theta[,index]
        sigma.q.C[[i]] = solve( 
          diag(1, Kp, Kp ) +
            ((A + J/2)/(b.q.sigma.me)) * (f_trace(Theta_i = Theta_i, Sig_q_Bpsi = sigma.q.Bpsi, Kp = Kp, Kt = Kt) + 
                                            t(mu.q.Bpsi) %*% Theta_i %*% t(Theta_i) %*% mu.q.Bpsi)
        )
        
        mu.q.C[i,] = ((A + J/2)/(b.q.sigma.me)) * sigma.q.C[[i]] %*% (psi.cur[,index]) %*% (Y[index] - fixef.est[index] )
      }
      
      for (i in 1:I) {
        index = which(subj.id == subjs[i])
        pcaef.est[index] = mu.q.C[i,] %*% psi.cur[,index]
      }
    }
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    if (Kp > 0) {
      ## measurement error variance
      resid = Y - fixef.est - pcaef.est
      b.q.sigma.me = as.numeric(B + .5 * (crossprod(resid) + 
                                            sum(diag(sumtXsXs %*% sigma.q.beta)) + 
                                            f_sum4(mu.q.c = mu.q.C, sig.q.c = sigma.q.C, mu.q.bpsi = mu.q.Bpsi, sig.q.bpsi = sigma.q.Bpsi, theta = Theta, subj.id = subj.id)) )
      
      ## lambda for FPCA basis functions
      for (K in 1:Kp) {
        b.q.lambda.Bpsi[K] = B + .5 * (t(mu.q.Bpsi[,K]) %*% diag(1, Kt, Kt) %*% mu.q.Bpsi[,K] + 
                                         sum(diag(diag(1, Kt, Kt) %*% sigma.q.Bpsi[(Kt*(K - 1) + 1):(Kt*K), (Kt*(K - 1) + 1):(Kt*K)])))
      }
    } else if (Kp == 0) {
      ## measurement error variance
      resid = Y - fixef.est - pcaef.est
      b.q.sigma.me = as.numeric(B + .5 * (crossprod(resid) + 
                                            sum(diag(sumtXsXs %*% sigma.q.beta)) ))
    }
    
    cat(".")
    
  }
  
  cat("\n")
  

  ## get coefficient functions over a common grid
  rownames(beta.cur) = c(trmstrings)
  beta.cur = t(beta.cur) %>% as.data.frame() %>%
    mutate(t = time) %>% 
    rename_(.dots = setNames(list("t"), time.var)) %>%
    arrange_(time.var) %>% unique() %>%
    subset(select = c(time.var, trmstrings))
  
  ## export fitted values
  Yhat = fixef.est + pcaef.est 
  
  ## export variance components
  sigeps.pm = 1 / as.numeric((A + J/2)/(b.q.sigma.me))
  
  ## export FPCA objects
  if (Kp > 0) {
    ## do svd to get rotated fpca basis
    rownames(psi.cur) = paste0("Psi.", 1:Kp)
    psi.cur = t(psi.cur) %>% as.data.frame() %>%
      mutate(t = time) %>% 
      arrange(t) %>% unique() %>%
      subset(select = c(paste0("Psi.", 1:Kp))) %>%
      as.matrix() %>%
      t()
    argvals = sort(unique(time))
    SVD = svd(psi.cur)
    efunctions = SVD$v
    scores = mu.q.C %*% SVD$u %*% diag(SVD$d, Kp, Kp)
    evalues = diag(cov(scores))
  } else if (Kp == 0) {
    argvals = sort(unique(time))
    efunctions = NULL
    scores = NULL
    evalues = NULL
  }
  
  Yhat.fpca = data.frame(index = time, value = pcaef.est, id = subj.id)
  class(Yhat.fpca) = c("data.frame", "refund_object")
  
  Y.fpca = data.frame(index = time, value = Y - fixef.est, id = subj.id)
  class(Yhat.fpca) = c("data.frame", "refund_object")
  
  fpca.obj = list(Yhat = Yhat.fpca,
                  Y = Y.fpca,
                  scores = scores,
                  mu = rep(0, length(argvals)),
                  efunctions = efunctions, 
                  evalues = evalues,
                  npc = Kp,
                  argvals = argvals)
  class(fpca.obj) = "fpca"
  
  ## export objects  
  ret = list(beta.cur, mu.q.beta, Yhat, sigeps.pm, fixef.est, 
             data.complete, data.model, formula, formula.model, mt_fixed, time.var, id.var, Kt, Kp, standardized,
             t.min, t.max, fpca.obj)
  names(ret) = c("beta.pm", "spline.coef.est", "Yhat", "sigeps.pm", "fixef",
                 "data", "data.model", "formula", "formula.model", "terms", "time.var", "id.var", "Kt", "Kp", "standardized",
                 "t.min", "t.max", "fpca.obj")
  
  class(ret) = "flcm"
  
  ret
  
}

