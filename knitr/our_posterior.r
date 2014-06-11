## This is a function that takes two parameters:
## - fitmodel: the model we want to evaluate
## - theta: the parameter value(s) at which we want to evaluate the posterior

our_posterior <- function(model, theta) {

    ## calculate the fitmodel prior for parameter vector theta using
    ## model$log.prior
    prior <- model$log.prior(theta)

    ## initialise the model (using model$initialise.theta and 
    ## the parameter vector theta)
    state.init <- model$initialise.state(theta)

    ## simulate the model at all the times given in model$data$time,
    ## that is at all the times that are in the data, using 
    ## model$simulate.model, with initial state given by staten.init
    trajectory <- model$simulate.model(theta, state.init, model$data$time)

    ## calculate the log-likelihood using model$log.likelihood
    likelihood <- model$log.likelihood(model$data, theta, trajectory)

    ## calulate the posterior using the log-prior and log-likelihood
    posterior <- prior + likelihood
    
    return(posterior)

}

