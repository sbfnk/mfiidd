Here is an example of how to code the posterior.


```r
## This is a function that takes two parameters:
## - fitmodel: the model we want to evaluate
## - theta: the parameter value(s) at which we want to evaluate the posterior

our_posterior <- function(model, theta) {

    ## INSERT HERE: calculate the fitmodel prior for parameter vector theta using
    ##              model$log.prior,and assign to a variable prior
    prior <- model$log.prior(theta)

    ## INSERT HERE: initialise the model (using model$initialise.theta and 
    ##              the parameter vector theta), and assign to a variable 
    ##              state.init
    state.init <- model$initialise.state(theta)

    ## INSERT HERE: simulate the model at all the times given in model$data$time,
    ##              that is at all the times that are in the data, using 
    ##              model$simulate.model, with initial state given by staten.init, 
    ##              and assign to a variable trajectory
    trajectory <- model$simulate.model(theta, state.init, model$data$time)

    ## INSERT HERE: calculate the log-likelihood using model$log.likelihood
    likelihood <- model$log.likelihood(model$data, theta, trajectory)

    ## INSERT HERE: calulate the posterior using the log-prior and log-likelihood
    posterior <- prior + likelihood
    
    return(posterior)

}
```

If you run into any problems, have a look at our [solution](mcmc_example_solution.md).
