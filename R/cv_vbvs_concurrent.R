#' cv_vbvs_concurrent
#' 
#' Cross validation to choose the tuning parameter v0 in the variational bayes variable selection algorithm
#' for the linear functional concurrent model. Uses five-fold cross validation.
#' 
#' @param Y matrix of functional responses
#' @param X list of matrices, each of which contains the functional predictor for an individual subject
#' @param Kt number of spline basis functions for coefficients and FPCs
#' @param Kp number of FPCs to estimate
#' @param v0 tuning parameter vector; normal spike variance. defaults to 0.01 to 0.1 in increments of 0.01.
#' @param v1 tuning parameter; normal slab variance
#' @param SEED seed value, used to ensure reproducibility of the training split
#' @param standardized logical; are covariates already standardized?
#' @param t.min minimum value to be evaluated on the time domain (useful if data are sparse and / or irregular). if `NULL`, taken to be minimum observed value.
#' @param t.max maximum value to be evaluated on the time domain (useful if data are sparse and / or irregular). if `NULL`, taken to be maximum observed value.
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' @export
cv_vbvs_concurrent = function(formula, id.var = NULL, data=NULL, Kt = 5, Kp = 2, v0 = seq(0.01, .1, .01), v1 = 100, 
                              SEED = 1, standardized = FALSE, t.min = NULL, t.max = NULL){
  
  OOS.sq.err = matrix(NA, nrow = length(v0), ncol = 5)
  outcome.var = as.character(formula)[2]
  set.seed(SEED)
  
  ## define folds
  subjs = unique(data[id.var][,1])
  groups = sample(subjs, length(subjs)) %>%
    split(., ceiling(seq_along(.)/ (length(.)/ 5)  ))
    
  for(FOLD in 1:5){

    data.train = filter(data, !(data[id.var][,1] %in% groups[[FOLD]]))
    data.test = filter(data, data[id.var][,1] %in% groups[[FOLD]])
  
    for(i.split in 1:length(v0)){
      cat(paste0(i.split, ". "))
      fit.vbvs = vbvs_concurrent(formula = formula, id.var = id.var, data = data.train, Kt = Kt, Kp = Kp, v0 = v0[i.split], v1 = v1, 
                                 standardized = standardized, t.min = t.min, t.max = t.max)
      fitted.test.vbvs = predict(fit.vbvs, data = data.test, standardized = standardized)
      OOS.sq.err[i.split, FOLD] = mean((data.test[outcome.var] - fitted.test.vbvs)^2, na.rm = TRUE)
    }
  }
  return(apply(OOS.sq.err, 1, mean))
  
}
