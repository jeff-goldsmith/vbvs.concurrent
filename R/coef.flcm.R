#' coef.flcm
#' 
#' Extract coefficient functions from a fitted concurrent functional linear model
#' 
#' This function is used to extract coefficient functions from a concurrent functional 
#' linear model, possibly on a user-specified grid.
#' 
#' @param object return object from \code{\link{vbvs_concurrent}} or
#'   \code{\link{vb_concurrent}}
#' @param t.new vector indicating the desired coordinates on which the functional 
#'   coefficients will be evaluated
#' @param ... these arguments are ignored
#'   
#' @return a data frame containing the evaluation points and coefficient function 
#'   values on the indicated grid
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' 
#' @import dplyr
#' @importFrom stats coef terms.formula quantile
#' @importFrom splines bs
#' 
#' @export
#' 
coef.flcm <- function(object, t.new = NULL, ...) {
  
  tf <- terms.formula(object$formula.model, specials = NULL)
  trmstrings <- attr(tf, "term.labels")

  time.var = object$time.var
  t.original = object$data.model[time.var][,1]
  
  if (is.null(t.new)) { t.new = t.original }
  
  if (min(t.new) < min(object$t.min) | max(t.new) > max(object$t.max)) {
    stop("Specified domain extends beyond domain on which model was originally estimated")
  }

  ## ensure that spline basis covers the same domain as the original; these points
  ## will be removed shortly
  t.new = c(object$t.min, object$t.max, t.new)
    
  Theta = t(bs(t.new, knots = quantile(t.original, probs = seq(0, 1, length = object$Kt - 2))[-c(1,object$Kt - 2)], 
               intercept = TRUE, degree = 3))
  beta.cur = t(object$spline.coef.est) %*% (Theta)
#  rownames(beta.cur) = c("int", trmstrings)
  rownames(beta.cur) = trmstrings
  beta.cur = t(beta.cur) %>% as.data.frame() %>%
    slice(., -(1:2)) %>%
    mutate(t = t.new[-(1:2)]) %>% 
#    subset(select= c("t", "int", trmstrings))
    select(t, one_of(trmstrings)) 
  
  names(beta.cur)[1] = time.var
  
  beta.cur
  
}
