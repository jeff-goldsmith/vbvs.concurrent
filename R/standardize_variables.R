#' standardize_variables
#' 
#' This function standardizes predictor and response variables prior to fitting the
#' concurrent model using variable selection.
#'
#' @param data.original data used to fit the model
#' @param data.new data to be standardized. by default, this is the same as 
#' \code{data.original}; however it may be different to allow predictions on new
#' datasets
#' @param time.var the variable name that indicates the time variable
#' @param trmstrings collection of variables to be standardized
#' @param LHS the response variable name
#' 
#' @importFrom FNN knn.index knnx.index
#' 
#' @export
#'
standardize_variables = function(data.original, data.new = NULL, time.var, trmstrings, LHS = NULL) {
  
  if (is.null(data.new)) { data.new = data.original}
  
  t.original = data.original[time.var][,1]
  t.new = data.new[time.var][,1]
  
  knn.index.original = as.data.frame(t(knnx.index(t.original, t.original, k = round(length(t.original)*.2))))
  knn.index.new = as.data.frame(t(knnx.index(t.original, t.new, k = round(length(t.original)*.2))))
  
  for (trm in trmstrings) {
    covar.cur = data.original[trm][,1]
    
    ## get mean and var for original data
    mean.fit.original = sapply(knn.index.original, function(u){mean(covar.cur[u])})
    sq.resid.original = (covar.cur - mean.fit.original) ^ 2
    
    ## standardize using computed mean and variance
    mean.fit.new = sapply(knn.index.new, function(u){mean(covar.cur[u])})
    var.fit.new = sapply(knn.index.new, function(u){mean(sq.resid.original[u])})
    
    data.new[trm] = (data.new[trm] - mean.fit.new) / sqrt(var.fit.new)
    cat(".")
  }
  
  for (trm in LHS) {
    covar.cur = data.original[trm][,1]
    
    ## get mean for original data
    mean.fit.new = sapply(knn.index.new, function(u){mean(covar.cur[u])})
    
    ## remove mean for new data
    data.new[trm] = (data.new[trm] - mean.fit.new)
    cat(".")
  }
  
  data.new
  
}