context("Prediction")

#### simulate example dataset for testing ####

library(dplyr)
library(tidyr)

## set design elements
set.seed(1)
I = 100
p = 50

## coefficient functions
beta1 = function(t) { sin(2*t*pi) }
beta2 = function(t) { cos(2*t*pi) }
beta3 = function(t) { 1 }

## FPC basis functions
psi1 = function(t) { sin(2*t*pi) }
psi2 = function(t) { cos(2*t*pi) }

## generate subjects, observation times, and FPC scores
time.data = sapply(1:I, function(u) {
  ji = sample(10:15, 1)
  rbind(runif(ji, 0, 1) %>% sort, 
        rep(u, ji),
        rep(rnorm(1, 0, 3), ji),
        rep(rnorm(1, 0, 1), ji))
}) %>% unlist() %>% matrix(ncol = 4, byrow = TRUE)
colnames(time.data) = c("time", "subj", "c_i1", "c_i2")
time.data = as.data.frame(time.data)

## generate predictor data
predictor.data = matrix(rnorm(dim(time.data)[1] * p), dim(time.data)[1], p)
colnames(predictor.data) = paste0("Cov_", 1:p)

## combine and generate responses
concurrent.data = cbind(time.data, predictor.data)
concurrent.data = 
  mutate(concurrent.data,
         Y = Cov_1 * beta1(time) +            ## fixed effects
           Cov_2 * beta2(time) + 
           Cov_3 * beta3(time) + 
           c_i1 * psi1(time) +              ## pca effects
           c_i2 * psi2(time) + 
           rnorm(dim(concurrent.data)[1]))  ## measurement error

pred.list = paste("Cov", 1:p, sep = "_")
formula = as.formula( paste("Y ~", paste(pred.list, collapse = "+"), "| time") )


#### tests ####

test_that("the predict function works for vb and vbvs", {
  
  fit.vbvs = vbvs_concurrent(formula, id.var = "subj", data = concurrent.data, standardized = TRUE,
                             t.min = 0, t.max = 1)
  fit.vb = vb_concurrent(formula, id.var = "subj", data = concurrent.data, standardized = TRUE,
                             t.min = 0, t.max = 1)
  
  expect_equal(max(abs(predict(fit.vbvs) - fit.vbvs$fixef)), 0, tolerance = 1e-5)
  expect_equal(max(abs(predict(fit.vb) - fit.vb$fixef)), 0, tolerance = 1e-5)
  
})

test_that("predict options work", {
  
  fit.vbvs = vbvs_concurrent(formula, id.var = "subj", data = concurrent.data, standardized = TRUE,
                             t.min = 0, t.max = 1)
  predictions = predict(fit.vbvs, data.new = filter(concurrent.data, subj == 3))
  
  expect_equal(max(abs(predictions - fit.vbvs$fixef[which(concurrent.data$subj == 3)])), 0, tolerance = 1e-5)

})
