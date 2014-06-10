#'Marginal log-likelihood for a deterministic model
#'
#'Compute the marginal log-likelihood of \code{theta} for a deterministic model defined in a \code{\link{fitmodel}} object.
#' @param theta named vector of model parameters. Names must correspond to those of \code{names(fitmodel$theta)}.
#' @param fitmodel \code{\link{fitmodel}} object.
#' @export
#' @return numeric value of the log-likelihood
marginalLogLikelihoodDeterministic <- function(theta, fitmodel) {

	data <- fitmodel$data

	# time sequence (must include initial time)
	times <- c(0,data$time)

	# simulate model at successive observation times of data
	traj <- fitmodel$simulate.model(theta,fitmodel$initialise.state(theta),times)

	# compute log-likelihood
	log.likelihood <- fitmodel$log.likelihood(data=data,theta=theta,model.traj=traj)

	return(log.likelihood)
}

#'Marginal log-likelihood for a stochastic model
#'
#'Compute a Monte-Carlo estimate of the log-likelihood of \code{theta} for a stochastic model defined in a \code{\link{fitmodel}} object, using \code{\link{bootstrapParticleFilter}}
#' @inheritParams marginalLogLikelihoodDeterministic
#' @inheritParams bootstrapParticleFilter
#' @export
#' @seealso bootstrapParticleFilter
#' @return Monte-Carlo estimate of the marginal log-likelihood of \code{theta}
marginalLogLikelihoodStochastic <- function(theta, fitmodel, n.particles, n.cores = 1) {

	# replace parameter values
	fitmodel$theta[names(theta)] <- theta

	# run SMC
	smc <- bootstrapParticleFilter(fitmodel=fitmodel, n.particles=n.particles, n.cores=n.cores)

	return(smc$log.likelihood)
}

#'Target posterior distribution for a fitmodel
#'
#'This function evaluates the posterior distribution at \code{theta} and returns the result in a suitable format for \code{\link{mcmcMH}}.
#' @param theta named vector of estimated theta
#' @param log.prior \R-function to compute the log-prior of \code{theta}, as returned by \code{\link{fitmodel}}
#' @param log.prior.args list of arguments passed to \code{log.prior}
#' @param marginal.log.likelihood \R-function to compute the marginal log-likelihood of \code{theta}
#' @param marginal.log.likelihood.args list of arguments passed to \code{marginal.log.likelihood}
#' @export
#' @return a list of two elements
#' \itemize{
#' 	\item \code{log.dist} numeric, logged value of the posterior distribution evaluated at \code{theta}
#' 	\item \code{trace} named vector with trace information (theta, log.prior, marginal.log.likelihood, log.posterior)
#' }
targetPosterior <- function(theta, log.prior, log.prior.args=list(), marginal.log.likelihood, marginal.log.likelihood.args=list()) {

	theta.log.prior <- do.call(log.prior, c(list(theta=theta),log.prior.args))

	if(is.finite(theta.log.prior)){
		theta.marginal.log.likelihood <- do.call(marginal.log.likelihood, c(list(theta=theta), marginal.log.likelihood.args))
	}else{
		# do not compute log.likelihood	(theta prior is 0)
		theta.marginal.log.likelihood  <-  -Inf
	}

	theta.log.posterior <- theta.log.prior + theta.marginal.log.likelihood

	return(list(log.dist=theta.log.posterior, trace=c(theta,log.prior=theta.log.prior,marginal.log.likelihood=theta.marginal.log.likelihood,log.posterior=theta.log.posterior)))

}



