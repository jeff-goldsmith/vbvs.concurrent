---
title: 'vbvs.concurrent: Fitting Methods for the Functional Linear Concurrent Model'
date: "16 December 2016"
tags:
- statistical analysis
- R
authors:
- affiliation: Columbia University
  name: Jeff Goldsmith
  orcid: 0000-0002-6150-8997
bibliography: paper.bib
---

# Summary

Functional data analysis is concerned with understanding measurements made over time, space, frequencies, and other domains for multiple subjects. Given the ubiquity of wearable devices, it is common to obtain several data streams monitoring blood pressure, physical activity, heart rate, location, and other quantities on study participants in parallel. Each of these data streams can be thought of as functional data, and the functional linear concurrent model is useful for relating predictor data to an outcome. This model can be written 
$$
Y_i(t) = \beta_0(t) + \sum_{k = 1}^{p}X_{ik}(t)\beta_k(t) + \delta_i(t)
$$
where $Y_i(t)$ is the functional response for subject $i$, the $X_{ik}(t)$ are functional predictors, the $\beta_k(t)$ are functional coefficients of interest, and the $\delta_i(t)$ are possibly correlated errors. This package implements two statistical methods (with and without variable selection) for estimating the parameters in the functional linear concurrent model; these methods are described in detail [@goldsmith2017]. 

Given tidy datasets containing functional responses and predictors for all subjects, `vb_concurrent` and `vbvs_concurrent` fit the functional linear concurrent model using variational Bayes and variational Bayes with variable selection, respectively. These functions produce objects of class `flcm` and have the same structure. Coefficients and predictions can be extracted or computed from `flmc` objects using `coef` and `predict`. Interactive visualizations of `flmc` objects are supported through the `refund.shiny` package [@refund.shiny, @wrobel2016]. 


# References
