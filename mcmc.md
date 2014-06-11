---
Title: Practical session: Fitting a deterministic model with MCMC
---

```
## Warning: no help found for 'fitmodel'
```

## Introduction

For this session, you will fit a deterministic SIR model to a recent measles outbreak in Europe, making use of the model simulation and likelihood functions from the last session. The data is in the "measles" data set, and can be loaded with


```r
data(measles)
```

You can look at the data set with


```r
measles
```

or visualise them using `plot` or `ggplot`, e.g.


```r
library(ggplot2)
ggplot(measles, aes(x = time, y = Inc)) + geom_point()
```

You will use the likelihood and model simulation functions from the last session.

## Prior

First, you need to decide on the prior distribution of the different parameters, such as `R0`, the infectious period, etc. The `fitparam` class has a structure set up for this: you can pass a prior distribution by setting the `prior` parameter. This expects an R `list` with two elements: `distribution`, a name of a density function for the prior distribution, and `parameter`, a named vector of parameters to pass to the `distribution` functions. Distribution functions in R start with the letter `d`. Common ones which might be useful for prior distributions are `dunif` and `dnorm`. A good resource for the most commonly used probability distributions can be found in the [R tutorial ebook](http://www.r-tutor.com/elementary-statistics/probability-distributions).

Example:

```r
R0 <- fitparam(name = "R0",
               value = 2,
               prior = list(distribution = "dunif",
                            parameters = c(min = 1, max = 5)))
```
To calculate the model prior, you have to take the product of all individual model prior (or sum the logged priors). There is already a function to do this for you, called `compositeLogPrior`. It takes a list of `fitparam` objects and returns the composite prior from the specified prior distributions.

Example: if you have fitparam objects `R0` and `infectious.period` defined,


```r
compositeLogPrior(list(R0))
```

```
## [1] -1.386
```

Using this, write a function to calculate the prior:

```r
my_logPrior <- function(list.fitparam) {

    ## calculate the composite prior for the given parameter list

    ## return the logged prior probability
}
```



Once you have done this, you can pass the prior function to a `fitmodel` object. For example, if you have two fitparam objects `R0` and `infectious.period`, and a function `my_logPrior`, you can call


```r
my_SIR <- fitmodel(name = "SIR",
                   state.variables = c("S", "I", "R"),
                   list.fitparam = list(R0),
                   log.prior.fitparam = my_logPrior)
```

```
## log.prior(theta) should return a single finite value
## Test: -1.386 
## --> log.prior looks good!
```

`fitmodel` converts the prior function into a function that takes a named parameter vector `theta`, available as `SIR$log.prior` (see the `fitmodel` documentation). You can then use this to evaluate the prior at different parameter values:


```r
my_SIR$log.prior(c(R0 = 3))
```

```
## [1] -1.386
```

You will probably get a different value, depending on your prior function.

## Posterior



To calculate the posterior for a given parameter, write a function that takes a vector `theta` of named parameters and returns a logged posterior probability.


```r
my_posterior <- function(fitmodel, theta) {

    ## calculate the fitmodel prior for parameter theta

    ## calculate the fitmodel likelihood for parameter theta

    ## return the logged posterior probability
}
```



When called with a `fitmodel` and a vector of named parameters `theta`, it should return the value of the posterior at that value of `theta`, e.g.


```r
theta <- c(R0 = 15, IP = 2, rho = 0.5, pI0 = 35 / 7.3e+6, pR0 = 0.9, N = 7.3e+6)
my_posterior(my_SIR, theta)
```

```
## [1] -103697
```

You will probably get a different value, depending on your prior and likelihood functions.

## MCMC

To sample from a distribution using MCMC, write a function that takes a target distribution like your posterior (i.e., a function that takes one argument, `theta`), an initial `init.theta`, a covariance matrix and the number of iterations. It should evaluate the target function at `init.theta` and then apply the Metropolis-Hastings algorithm for the specified number of iterations.

To draw a random vector from a multivariate Gaussian proposal distribution, you can use the function `rmvnorm`, see `?rmvnorm`. To draw a random number between 0 and 1, you can use `runif(n = 1)`.


```r
my_mcmc <- function(target.dist, init.theta, covariance.matrix, n.iterations) {

    ## evaluate the function target_dist at init_theta

    ## repeat n_iterations times:

    ## - draw a new theta from the (multivariate Gaussian) proposal distribution

    ## - evaluate the function target_dist at the proposed theta

    ## - draw a random number between 0 and 1

    ## - accept or reject by comparing the random number to the acceptance probability

    ## return the chain (i.e., the value of the current theta at every iteration)

}
```

You will probably find it useful to save the state of the chain (i.e., the value of `theta` at every iteration) and the current acceptance probability (i.e., the proportion of proposed `theta`s that are being accepted) in a data frame or vector.

To get the covariance matrix, you can use the `sd.proposal` argument of `fitparam`, e.g. This sets the covariance matrix to be diagonal, that is the proposals for each parameter will not depend on the other parameters, but will be independently normally distributed, with the standard deviation given by `sd.proposal`.


```r
R0 <- fitparam(name = "R0",
               value = 2,
               prior = list(distribution = "dunif",
                            parameters = c(min = 1, max = 5)),
               sd.proposal = 0.1)
```

Once you pass this to a `fitmodel` (in `parameters`), it will return the covariance matrix containing all the parameters in `fitmodel$gaussian.proposal$covmat`.




```r
my_SIR <- fitmodel(name = "SIR",
                   state.variables = c("S", "I", "R", "Inc"),
                   list.fitparam = list(R0, InfectiousPeriod, ReportingRate, PopSize,
                                        proportionI0, proportionR0),
                   log.prior.fitparam = my_logPrior,
                   log.likelihood = my_logLikelihood,
                   simulate.model = my_simulateDeterministic,
                   initialise.state = my_initialiseState,
                   data = measles,
                   verbose = F)
my_SIR$gaussian.proposal$covmat
```

```
##       R0   IP  rho N       pI0   pR0
## R0  0.01 0.00 0.00 0 0.000e+00 0e+00
## IP  0.00 0.25 0.00 0 0.000e+00 0e+00
## rho 0.00 0.00 0.01 0 0.000e+00 0e+00
## N   0.00 0.00 0.00 0 0.000e+00 0e+00
## pI0 0.00 0.00 0.00 0 1.877e-14 0e+00
## pR0 0.00 0.00 0.00 0 0.000e+00 1e-04
```

You will probably see a different matrix, depending on the proposal standard deviations you have set.
