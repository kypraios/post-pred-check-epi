# post-pred-check-epi

Posterior Predictive Checking for Partially Observed Stochastic Epidemic Models

## Overview

This repository contains `R` code for implementing the methods presented in

>Aristotelous G., Kypraios, T. and O'Neil, P.D. (2022) Posterior Predictive Checking for Partially Observed Stochastic Epidemic Models, *Bayesian Analysis*; to appear.  

which addresses the problem of assessing the fit of stochastic epidemic models to data. Two novel model assessment methods are developed and presented in the manuscript above, based on disease progression curves, namely the **distance method** and the **position-time method**. Both methods provide visual and quantitative outputs with meaningful interpretations. 

We illustrate both methods by applying them to data from an outbreak which occurred in a boarding school in England and led to 102 cases among 1307 students [Smith *et al.*, 2009](https://www.eurosurveillance.org/content/10.2807/ese.14.27.19263-en).

##  Usage 

The file `functions.R` contains several `R` functions, two of which are the `distance.method` and `position.time.method` and they correspond to the two different methods for assessing the goodness of fit of a stochastic epidemic model when fitted to data. Both function require as arguments a vector with the observed removal curve (i.e. recovery times) and a matrix containing samples from the posterior predictive distribution of the removal times once the model has been fitted. 

The two methods can be implemented by calling these two functions as follows:

`distance.method(R.obs, R.STAR)`

`position.time.method(R.obs, R.STAR)`

where `R.obs` denotes the vector with the observed removal curves and `R.STAR` denotes the sampled removal curves.

We reproduce some the results in Section 8 in the `illustration.Rmd`
