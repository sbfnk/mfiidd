Below, you can find an example of how to code the posterior. Some bits are left out for you to fill in (marked "INSERT HERE"). Each "INSERT HERE" statement requires one line of code. If you struggle, you can find a link to the solution below the function.


```r
## This is a function that takes two parameters:
## - model: the model we want to evaluate
## - theta: the parameter value(s) at which we want to evaluate the posterior

our_posterior <- function(model, theta) {

    ## INSERT HERE: calculate the fitmodel prior for parameter vector theta using
    ##              model$log.prior,and assign to a variable prior

    ## INSERT HERE: initialise the model (using model$initialise.theta and 
    ##              the parameter vector theta), and assign to a variable 
    ##              state.init

    ## INSERT HERE: simulate the model at all the times given in model$data$time,
    ##              that is at all the times that are in the data, using 
    ##              model$simulate.model, with initial state given by staten.init, 
    ##              and assign to a variable trajectory

    ## INSERT HERE: calculate the log-likelihood using model$log.likelihood

    ## INSERT HERE: calulate the posterior using the log-prior and log-likelihood
    
    return(posterior)

}
```

If you run into any problems, have a look at our [solution](posterior_example_solution.md).
