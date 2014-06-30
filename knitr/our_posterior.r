## @knitr posterior
## This is a function that takes 4 arguments:
## - fitmodel, a fitmodel object that defines the model dynamics,
##   prior and likelihoods.
## - theta, a named vector of parameters
## - state.init,  a named vector of initial conditions
## - data, the data set we are fitting the model to
## It should return the posterior for the given model, parameters,
## initial conditions and data.
our_posterior <- function(fitmodel, theta, state.init, data) {

    ## calculate the fitmodel prior for parameter vector theta using
    ## fitmodel$logPrior, and assign to variable log.prior
    log.prior <- model$logPrior(theta)

    ## calculate the fitmodel prior for parameter vector theta using
    ## fitmodel$logPrior, and assign to variable log.prior
    log.likelihood <- trajLogLike(fitmodel, theta, state.init, data)

    ## calulate the posterior using the log-prior and log-likelihood
    posterior <- log.prior + log.likelihood

    return(posterior)

}
