#'Compute the distance to the data for ABC
#'
#'Compute the distance (using \code{fitmodel$distance.ABC}) between the observed time series and a simulated time series obtained by running the model with parameters \code{theta}.
#' @inheritParams marginalLogLikelihoodDeterministic
#' @export
#' @return numeric value of the log-likelihood
computeDistanceABC <- function(theta, fitmodel) {

	data <- fitmodel$data

	# time sequence (must include initial time)
	times <- c(0,data$time)

	# simulate model at successive observation times of data
	traj <- fitmodel$simulate.model(theta,fitmodel$initialise.state(theta),times)

	# generate simulated observation
	traj.obs <- fitmodel$generate.observation(model.traj=traj,theta=theta)

	# compute distance
	dist.ABC <- fitmodel$distance.ABC(data=data,model.traj.obs=traj.obs)

	return(dist.ABC)
}

#'Target ABC posterior distribution for a fitmodel
#'
#'This function evaluates the ABC posterior distribution at \code{theta} and returns the result in a suitable format for \code{\link{mcmcMH}}.
#' @param epsilon numeric vector, ABC tolerances for distances between data and simulations.
#' @inheritParams targetPosterior
#' @inheritParams marginalLogLikelihoodDeterministic
#' @export
#' @seealso computeDistanceABC
#' @return a list of two elements
#' \itemize{
#' 	\item \code{log.dist} numeric, logged value of the ABC posterior distribution evaluated at \code{theta}
#' 	\item \code{trace} named vector with trace information (theta, log.prior, distance.ABC, log.posterior)
#' }
targetPosteriorABC <- function(theta,fitmodel,epsilon) {

	theta.log.prior <- fitmodel$log.prior(theta=theta)

	if(is.finite(theta.log.prior)){
		theta.dist.ABC <- computeDistanceABC(theta,fitmodel)
		if(length(epsilon)!=length(theta.dist.ABC)){
			stop("Length of ",sQuote("epsilon")," differs from the length of the distance vector returned by ",sQuote("fitmodel$distance.ABC"),call.=FALSE)
		}
		theta.log.posterior <- theta.log.prior + log(ifelse(all(theta.dist.ABC <= epsilon),1,0))

	}else{
		# do not compute ABC distance (theta prior is 0)
		theta.dist.ABC  <- Inf
		theta.log.posterior  <-  -Inf
	}


	return(list(log.dist=theta.log.posterior, trace=c(theta,log.prior=theta.log.prior,distance.ABC=theta.dist.ABC,log.posterior=theta.log.posterior)))

}
