#' predict.flcm
#' 
#' Compute predicted values for a new dataset based on a fitted concurrent functional linear model
#' 
#' @param object return object from \code{\link{vbvs_concurrent}} or
#'   \code{\link{vb_concurrent}}
#' @param data.new dataset for which fitted values should be computed. should be a dataframe with names the same
#'   as those used to fit the model, including the time variable used to parameterize functional observations
#' @param standardized logical; are values in `data.new` standardized? if not, they will be standardized using the
#'   same scheme as used in the original analysis
#' @param ... these arguments are ignored
#' 
#' @return a vector of fitted values
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' 
#' @importFrom stats predict terms.formula
#' 
#' @export
#' 
predict.flcm <- function(object, data.new = NULL, standardized = TRUE, ...) {
  
  if (is.null(data.new)) {data.new = object$data.model}
  
  ## add check that column names match?
  
  tf <- terms.formula(object$formula.model, specials = NULL)
  trmstrings <- attr(tf, "term.labels")
  
  ## normalize variables
  
  if (!standardized) {
    
    cat("Standardizing Variables \n")
    data.new = standardize_variables(data.original = object$data, data.new = data.new, 
                                       time.var = object$time.var, trmstrings = trmstrings, LHS = NULL)
    cat("\n")
    
  }
  
  ## get coefficients for time on new dataset
  beta.hat = coef(object, t.new = data.new[object$time.var][,1])
  
  ## multiply each coefficient by each observed value; sum to get fitted value
#  fitted = mutate(data.new, int = 1) %>% subset(select = c("int", trmstrings)) %>% 
  fitted = data.new %>% subset(select = trmstrings) %>% 
    '*'(., select(beta.hat, -contains(object$time.var))) %>%
    apply(1, sum)

  fitted 
  
}
