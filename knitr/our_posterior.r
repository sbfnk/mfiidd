## @knitr posterior
# This is a function that takes 4 arguments:
# - fitmodel, a fitmodel object that defines the model dynamics,
#   prior and likelihoods.
# - theta, a named vector of parameters
# - init.state,  a named vector of initial conditions
# - data, the data set we are fitting the model to
# It returns the posterior for the given model, parameters, initial
# conditions and data.
my_logPosterior <- function(fitmodel, theta, init.state, data) {

    # calculate the fitmodel prior for parameter vector theta using
    # fitmodel$logPrior, and assign to variable log.prior
    log.prior <- fitmodel$logPrior(theta)

    # calculate the fitmodel prior for parameter vector theta using
    # fitmodel$logPrior, and assign to variable log.prior
    log.likelihood <- trajLogLike(fitmodel, theta, init.state, data)

    # calulate the log-posterior using the log-prior and log-likelihood
    log.posterior <- log.prior + log.likelihood

    return(log.posterior)

}
