#' cv_vbvs_concurrent
#' 
#' Cross validation to choose the tuning parameter v0 in the variational bayes variable selection algorithm
#' for the linear functional concurrent model. Uses five-fold cross validation.
#' 
#' @param formula formula for desired regression. should have form \code{ y ~ x1 + x2 + ... + x_k | t}, where \code{t} is the variable that parameterized observed functions
#' @param id.var variable giving subject ID vector
#' @param data optional data frame
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
#' 
#' @importFrom FNN knn.index knnx.index
#' 
#' @export
cv_vbvs_concurrent = function(formula, id.var = NULL, data = NULL, Kt = 5, Kp = 2, v0 = seq(0.01, .1, .01), v1 = 100, 
                              SEED = 1, standardized = FALSE, t.min = NULL, t.max = NULL){
  
  OOS.sq.err = matrix(NA, nrow = length(v0), ncol = 5)
  outcome.var = as.character(formula)[2]
  set.seed(SEED)
  
  ## define folds
  subjs = unique(pull(data, id.var))
  groups = sample(subjs, length(subjs)) %>%
    split(., ceiling(seq_along(.) / (length(.) / 5)  ))
    
  for (FOLD in 1:5) {

    data.train = filter(data, !(pull(data, id.var) %in% groups[[FOLD]]))
    data.test = filter(data, pull(data, id.var) %in% groups[[FOLD]])
  
    for (i.split in 1:length(v0)) {
      cat(paste0(i.split, ". "))
      fit.vbvs = vbvs_concurrent(formula = formula, id.var = id.var, data = data.train, Kt = Kt, Kp = Kp, v0 = v0[i.split], v1 = v1, 
                                 standardized = standardized, t.min = t.min, t.max = t.max)
      fitted.test.vbvs = predict(fit.vbvs, data = data.test, standardized = standardized)
      
      ## if variables aren't standardised, need to standardize response in data.test using values in data.train
      
      if (!standardized) {
        time.var = fit.vbvs$time.var
        
        t.original = pull(data.train, time.var)
        t.new = pull(data.test, time.var)
        knn.index.new = as.data.frame(t(knnx.index(t.original, t.new, k = round(length(t.original)*.2))))
        
        covar.cur = pull(data.train, outcome.var)
        
        mean.fit.new = sapply(knn.index.new, function(u){mean(covar.cur[u])})
        
        data.test[outcome.var] = (data.test[outcome.var] - mean.fit.new)
      }
      
      OOS.sq.err[i.split, FOLD] = mean((data.test[outcome.var] - fitted.test.vbvs) ^ 2, na.rm = TRUE)
    }
  }
  return(apply(OOS.sq.err, 1, mean))
  
}
