#' predict.concurFLM
#' 
#' Compute predicted values for a new dataset based on a fitted concurrent functional linear model
#' 
#' @param object return object from \code{\link{vbvs_concurrent}}, 
#'   \code{\link{vb_concurrent}}, or \code{\link{ols_concurrent}}
#' @param data.new dataset for which fitted values should be computed. should be a dataframe with names the same
#'   as those used to fit the model, including the time variable used to parameterize functional observations
#' @param standardized logical; are values in `data.new` standardized? if not, they will be standardized using the
#'   same scheme as used in the original analysis
#'   
#' @return a vector of fitted values
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' 
#' @importFrom stats predict
#' @importFrom FNN knn.index knnx.index
#' 
#' @export
#' 
predict.concurFLM <- function(object, data.new = NULL, standardized = FALSE) {
  
  ## add check that column names match?
  
  tf <- terms.formula(object$formula.model, specials = NULL)
  trmstrings <- attr(tf, "term.labels")
  
  ## normalize variables
  
  if(!standardized){
    cat("Standardizing Variables \n")
    
    t.original = object$data[object$time.var][,1]
    t.new = data.new[object$time.var][,1]

    knn.index.original = as.data.frame(t(knnx.index(t.original, t.original, k = round(length(t.original)*.2))))
    knn.index.new = as.data.frame(t(knnx.index(t.original, t.new, k = round(length(t.original)*.2))))
    
    for(trm in trmstrings){
      covar.cur = object$data[trm][,1]
      
      ## get mean and var for original data
      mean.fit.original = sapply(knn.index.original, function(u){mean(covar.cur[u])})
      sq.resid.original = (covar.cur - mean.fit.original)^2

      ## standardize using computed mean and variance
      mean.fit.new = sapply(knn.index.new, function(u){mean(covar.cur[u])})
      var.fit.new = sapply(knn.index.new, function(u){mean(sq.resid.original[u])})
      
      data.new[trm] = (data.new[trm] - mean.fit.new) / sqrt(var.fit.new)
      cat(".")
    }
    cat("\n")
  }
  
  ## get coefficients for time on new dataset
  beta.hat = coef(object, t.new = data.new[object$time.var][,1])
  
  ## multiply each coefficient by each observed value; sum to get fitted value
  fitted = mutate(data.new, int = 1) %>% subset(select = c("int", trmstrings)) %>% 
    '*'(., subset(beta.hat, select = -c(t))) %>%
    apply(1, sum)
  
  fitted 
  
}
