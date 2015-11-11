#' cv_vbvs_concurrent
#' 
#' Cross validation to choose the tuning parameter v0 in the variational bayes variable selection algorithm
#' for the linear functional concurrent model
#' 
#' @param Y matrix of functional responses
#' @param X list of matrices, each of which contains the functional predictor for an individual subject
#' @param Kt number of spline basis functions for coefficients and FPCs
#' @param Kp number of FPCs to estimate
#' @param v0 tuning parameter vector; normal spike variance. defaults to 0.01 to 0.1 in increments of 0.01.
#' @param v1 tuning parameter; normal slab variance
#' @param train.prop proportion of subjects to include in the training set
#' @param SEED seed value, used to ensure reproducibility of the training split
#' @param standardized logical; are covariates already standardized?
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). if `NULL`, taken to be minium observed value.
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). if `NULL`, taken to be minium observed value.
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' @export
cv_vbvs_concurrent = function(formula, id.var = NULL, data=NULL, Kt = 5, Kp = 2, v0 = seq(0.01, .1, .01), v1 = 100, 
                              train.prop = .8, SEED = 1, standardized = FALSE, t.min = NULL, t.max = NULL){
  
  OOS.sq.err = rep(NA, nrow = length(v0))
  outcome.var = as.character(form)[2]
  set.seed(SEED)
  
  ## split data
  subjs = unique(data[id.var][,1])
  in.sample = sample(subjs, length(subjs) * train.prop)
  
  data.train = filter(data, data[id.var][,1] %in% in.sample)
  data.test = filter(data, !(data[id.var][,1] %in% in.sample))
  
  for(i.split in 1:length(v0)){
    cat(paste0(i.split, ". "))
    fit.vbvs = vbvs_concurrent(formula = formula, id.var = id.var, data = data.train, Kt = Kt, Kp = Kp, v0 = v0[i.split], v1 = v1, 
                               standardized = standardized, t.min = t.min, t.max = t.max)
    fitted.test.vbvs = predict(fit.vbvs, data = data.test, standardized = standardized)
    OOS.sq.err[i.split] = mean((data.test[outcome.var] - fitted.test.vbvs)^2, na.rm = TRUE)
  }

  return(OOS.sq.err)
  
}