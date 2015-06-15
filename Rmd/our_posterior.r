## @knitr posterior
# This is a function that takes 4 arguments:
# - fitmodel, a fitmodel object that defines the model dynamics,
#   prior and likelihoods.
# - theta, a named vector of parameters
# - init.state,  a named vector of initial state
# - data, the data set we are fitting the model to
# It returns the posterior for the given model, parameters, initial
# state and data.
my_dLogPosterior <- function(fitmodel, theta, init.state, data) {

    # calculate the fitmodel prior for parameter vector theta using
    # fitmodel$dprior, and assign to variable log.prior
	log.prior <- fitmodel$dprior(theta, log = TRUE)

	# calculate the log-likelihood of `theta`
    # and `init.state` with respect to the data using `dTrajObs`
    # and assign to a variable `log.likelihood`    
    log.likelihood <- dTrajObs(fitmodel, theta, init.state, data, log = TRUE)

    # calulate the log-posterior using the log-prior and log-likelihood
	log.posterior <- log.prior + log.likelihood

	return(log.posterior)

}
