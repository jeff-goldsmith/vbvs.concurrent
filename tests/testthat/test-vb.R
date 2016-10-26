context("VB")

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

test_that("vb works on toy example", {
  
  fit.vb = vb_concurrent(formula, id.var = "subj", data = concurrent.data, standardized = TRUE,
                             t.min = 0, t.max = 1)
  
  expect_equal(dim(fit.vb$beta.pm)[1], dim(concurrent.data)[1])
  expect_equal(dim(fit.vb$beta.pm)[2], p + 1)
  expect_equal(mean((fit.vb$beta.pm$Cov_1 - beta1(fit.vb$beta.pm$time)) ^ 2), 0, tol = .1)
  
})

test_that("vb options work", {
  
  fit.vb = vb_concurrent(formula, id.var = "subj", data = concurrent.data, standardized = TRUE,
                             t.min = 0, t.max = 1, Kp = 1, Kt = 10)
  
  expect_equal(dim(fit.vb$spline.coef.est)[1], 10)
  expect_equal(dim(fit.vb$fpca.obj$scores)[2], 1)
  
  fit.vb = vb_concurrent(formula, id.var = "subj", data = concurrent.data, standardized = TRUE,
                             t.min = 0, t.max = 1, v1 = 0.0001)
  
  expect_equal(mean((fit.vb$beta.pm$Cov_1 - 0) ^ 2), 0, tol = .1)
  
})

test_that("vb returns objects of appropriate class", {
  
  fit.vb = vb_concurrent(formula, id.var = "subj", data = concurrent.data, standardized = TRUE,
                             t.min = 0, t.max = 1, Kp = 1, Kt = 10)
  
  expect_is(fit.vb, "flcm")
  expect_is(fit.vb$fpca.obj, "fpca")

})