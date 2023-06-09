## @knitr posterior
# This is a function that takes 4 arguments:
# - fitmodel, a fitmodel object that defines the model dynamics,
#   prior and likelihoods.
# - theta, a named vector of parameters
# - initState,  a named vector of initial state
# - data, the data set we are fitting the model to
# It returns the posterior for the given model, parameters, initial
# state and data.
my_dLogPosterior <- function(fitmodel, theta, initState, data) {
  # calculate the fitmodel prior for parameter vector theta using
  # fitmodel$dPrior, and assign to variable logPrior
  logPrior <- fitmodel$dPrior(theta, log = TRUE)

  # calculate the log-likelihood of `theta`
  # and `initState` with respect to the data using `dTrajObs`
  # and assign to a variable `logLikelihood`
  logLikelihood <- dTrajObs(fitmodel, theta, initState, data, log = TRUE)

  # calulate the log-posterior using the log-prior and log-likelihood
  logPosterior <- logPrior + logLikelihood

  return(logPosterior)
}
