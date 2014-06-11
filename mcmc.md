# Practical session: Fitting a deterministic model with MCMC



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
## This is a function that takes four parameters:
## - fitmodel: the model we want to evaluate
## - theta: the parameter value(s) at which we want to evaluate the posterior
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
## Error: object 'likelihod' not found
```

You will probably get a different value, depending on your prior and likelihood functions.

If you have trouble filling the empty bits in, have a look at our [example](posterior_example.md).

## MCMC

To sample from a distribution using MCMC, write a function that takes a target distribution like your posterior (i.e., a function that takes one argument, `theta`), an initial `init.theta`, a covariance matrix and the number of iterations. It should evaluate the target function at `init.theta` and then apply the Metropolis-Hastings algorithm for the specified number of iterations.

To draw a random vector from a multivariate Gaussian proposal distribution, you can use the function `rmvnorm`, see `?rmvnorm`. To draw a random number between 0 and 1, you can use `runif(n = 1)`.


```r
## This is a function that takes four parameters:
## - target: the target distribution, a function that takes one argument
##           (a vector) and returns the (logged) value of a distribution
## - theta.init: the initial value of theta, a named vector
## - covmat.proposal: the covariance matrix of the (Gaussian) proposal distribution,
##                    in the same order as in the "target" vector
## - n.iterations: the number of iterations
## it returns an MCMC trace (value of theta and target(theta) at every MCMC step)
my_mcmc <- function(target, theta.init, covmat.proposal, n.iterations) {

    ## evaluate the function "target" at "theta.init"

    ## repeat n.iterations times:

    ## - draw a new theta from the (multivariate Gaussian) proposal distribution

    ## - evaluate the function target at the proposed theta

    ## - draw a random number between 0 and 1

    ## - accept or reject by comparing the random number to the acceptance probability
    ##   use a Gaussian proposal distribution; what does this look like for the 
    ##   multivariate Gaussian? It's easiest if you assume the target distribution
    ##   returns the logarithm of its value at theta.

    ## return the chain (i.e., the value of the current theta at every iteration)

}
```



You will probably find it useful to save the state of the chain (i.e., the value of `theta` at every iteration) and the current acceptance probability (i.e., the proportion of proposed `theta`s that are being accepted) in a data frame or vector.

If you have trouble filling the empty bits in, have a look at our [example](mcmc_example.md).

You can use the Metropolis-Hastings sampler to sample from any target distribution. For example, imagine you don't know how to draw a normal distribution. You could use your sampler for this. First, let's define an intermediate function that returns the log of the N(0, 1) normal probability distribution at a given sampled value:


```r
dnorm.log <- function(x) {
  return(dnorm(x, log = T))
}
```

Next, use your MCMC sampler

```r
sample.sd <- 1
sigma <- matrix(sample.sd, nrow = 1)

starting.value <- c(x = 1)
iter <- 1000
trace <- my_mcmc(target = dnorm.log, theta.init = starting.value, 
                 covmat.proposal = sigma, n.iterations = iter)
```

You can visualise the sampling result using


```r
library(ggplot2)
p <- ggplot(trace, aes(x = x, y = ..count../sum(..count..)))
p <- p + geom_histogram(binwidth = 0.2)
p
```

and a trace of the sampler using


```r
trace$iteration <- seq(1, nrow(trace))

p <- ggplot(trace, aes(x = iteration, y = x))
p <- p + geom_line()
p
```

## Using MCMC to sample from the posterior

To get the covariance matrix of a `fitmodel`, you can use the `sd.proposal` argument of `fitparam`, e.g. This sets the covariance matrix to be diagonal, that is the proposals for each parameter will not depend on the other parameters, but will be independently normally distributed, with the standard deviation given by `sd.proposal`.


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

Once you have this, you can pass the proposal covariance matrix and posterior distributions to your MCMC sampler to explore the posterior. For this, create a wrapper function for the posterior:


```r
my_SIR_posterior <- function(theta) {
    my_posterior(my_SIR, theta)
}
```

and run MCMC on that function.

## Further topics

Have a look at the function `mcmcMH` by typing

```r
mcmcMH
```

```
## function (target, target.args, theta.init, gaussian.proposal = list(covmat = NULL, 
##     lower = NULL, upper = NULL), n.iterations, adapt.size.start = n.iterations, 
##     adapt.size.cooling = 0.99, adapt.shape.start = n.iterations, 
##     print.info.every = n.iterations/100) 
## {
##     theta.current <- theta.init
##     theta.propose <- theta.init
##     covmat.proposal <- gaussian.proposal$covmat
##     lower.proposal <- gaussian.proposal$lower
##     upper.proposal <- gaussian.proposal$upper
##     theta.names <- names(theta.init)
##     covmat.proposal <- covmat.proposal[theta.names, theta.names]
##     lower.proposal <- lower.proposal[theta.names]
##     upper.proposal <- upper.proposal[theta.names]
##     covmat.proposal.init <- covmat.proposal
##     start.adapt.size <- TRUE
##     start.adapt.shape <- TRUE
##     theta.estimated.names <- names(which(diag(covmat.proposal) > 
##         0))
##     target.theta.current <- do.call(target, c(list(theta = theta.current), 
##         target.args))
##     trace <- data.frame(t(target.theta.current$trace), weight = 1)
##     acceptance.rate <- 0
##     scaling.sd <- 1
##     covmat.empirical <- covmat.proposal
##     covmat.empirical[, ] <- 0
##     theta.mean <- theta.current
##     if (is.null(print.info.every)) {
##         print.info.every <- n.iterations + 1
##     }
##     start_iteration_time <- Sys.time()
##     for (i.iteration in seq_len(n.iterations)) {
##         if (i.iteration >= adapt.size.start && acceptance.rate * 
##             i.iteration < adapt.shape.start) {
##             if (start.adapt.size) {
##                 message("\n---> Start adapting size of covariance matrix")
##                 start.adapt.size <- 0
##             }
##             scaling.sd <- scaling.sd * exp(adapt.size.cooling^(i.iteration - 
##                 adapt.size.start) * (acceptance.rate - 0.234))
##             covmat.proposal <- scaling.sd^2 * covmat.proposal.init
##         }
##         else if (acceptance.rate * i.iteration >= adapt.shape.start) {
##             if (start.adapt.shape) {
##                 message("\n---> Start adapting shape of covariance matrix")
##                 start.adapt.shape <- 0
##             }
##             covmat.proposal <- 2.38^2/length(theta.estimated.names) * 
##                 covmat.empirical
##         }
##         if (i.iteration%%round(print.info.every) == 0) {
##             end_iteration_time <- Sys.time()
##             x <- trace[nrow(trace), ]
##             x <- paste(paste0(names(x), "=", sprintf("%.2f", 
##                 x)), collapse = "|")
##             suppressMessages(time.estimation <- round(as.period((end_iteration_time - 
##                 start_iteration_time) * 10000/round(print.info.every))))
##             message("Iteration: ", i.iteration, "/", n.iterations, 
##                 " Time 10000 iter: ", time.estimation, " Acceptance rate: ", 
##                 sprintf("%.3f", acceptance.rate), " Scaling.sd: ", 
##                 sprintf("%.3f", scaling.sd), " State:", x)
##             start_iteration_time <- end_iteration_time
##         }
##         if (any(diag(covmat.proposal)[theta.estimated.names] < 
##             .Machine$double.eps)) {
##             print(covmat.proposal[theta.estimated.names, theta.estimated.names])
##             stop("non-positive definite covmat", call. = FALSE)
##         }
##         theta.propose[theta.estimated.names] <- as.vector(rtmvnorm(1, 
##             mean = theta.current[theta.estimated.names], sigma = covmat.proposal[theta.estimated.names, 
##                 theta.estimated.names], lower = lower.proposal[theta.estimated.names], 
##             upper = upper.proposal[theta.estimated.names]))
##         target.theta.propose <- do.call(target, c(list(theta = theta.propose), 
##             target.args))
##         if (!is.finite(target.theta.propose$log.dist)) {
##             log.acceptance <- -Inf
##         }
##         else {
##             log.acceptance <- target.theta.propose$log.dist - 
##                 target.theta.current$log.dist + dtmvnorm(x = theta.current[theta.estimated.names], 
##                 mean = theta.propose[theta.estimated.names], 
##                 sigma = covmat.proposal[theta.estimated.names, 
##                   theta.estimated.names], lower = lower.proposal[theta.estimated.names], 
##                 upper = upper.proposal[theta.estimated.names], 
##                 log = TRUE) - dtmvnorm(x = theta.propose[theta.estimated.names], 
##                 mean = theta.current[theta.estimated.names], 
##                 sigma = covmat.proposal[theta.estimated.names, 
##                   theta.estimated.names], lower = lower.proposal[theta.estimated.names], 
##                 upper = upper.proposal[theta.estimated.names], 
##                 log = TRUE)
##         }
##         if (is.accepted <- (log(runif(1)) < log.acceptance)) {
##             trace <- rbind(trace, c(target.theta.propose$trace, 
##                 weight = 1))
##             theta.current <- theta.propose
##             target.theta.current <- target.theta.propose
##         }
##         else {
##             trace$weight[nrow(trace)] <- trace$weight[nrow(trace)] + 
##                 1
##         }
##         acceptance.rate <- acceptance.rate + (is.accepted - acceptance.rate)/i.iteration
##         tmp <- updateCovmat(covmat.empirical, theta.mean, theta.current, 
##             i.iteration)
##         covmat.empirical <- tmp$covmat
##         theta.mean <- tmp$theta.mean
##     }
##     return(list(trace = trace, acceptance.rate = acceptance.rate, 
##         covmat.empirical = covmat.empirical))
## }
## <environment: namespace:fitcourseR>
```

Can you understand what the function does? You will see that the function uses `rmvtnorm` and `dmvtnorm`, truncated Gaussian distributions, for the proposal. Can you think of why it does that?
