#'Metropolis-Hasting MCMC
#'
#'Run \code{n.iterations} of a Metropolis-Hasting MCMC to sample from the target distribution using a gaussian proposal kernel, possibly truncated.
#' @param target \R-function that takes at least one argument, called \code{theta}, which is the named vector of parameter of the model. This function returns a list of 2 elements:
#' \itemize{
#' \item \code{log.dist} the logged value of the target distribution, evaluated at \code{theta}.
#' \item \code{trace} a named vector of value that will be printed in the \code{trace} data.frame returned by \code{mcmcMH}.
#' }
#' @param target.args list of additional arguments to be passed to the function \code{target}.
#' @param theta.init named vector of initial parameter values to start the chain.
#' @param gaussian.proposal list of theta for the - potentially truncated - multi-variate normal proposal kernel of the MCMC. Contains 3 elements:
#'\itemize{
#'	\item \code{covmat} covariance matrix. Must have named rows and columns with at least all estimated theta.
#'	\item \code{lower} named vector of lower truncation points.
#'	\item \code{upper} named vector of upper truncation points.
#'} 
#' @param n.iterations number of iterations to run the MCMC chain.
#' @param adapt.size.start number of iterations to run before adapting the size of the proposal covariance matrix (see note below).
#' @param adapt.size.cooling cooling factor for the scaling factor of the covariance matrix during size adaptation (see note below).
#' @param adapt.shape.start number of accepted jumps before adapting the shape of the proposal covariance matrix (see note below).
#' @param print.info.every frequency of information on the chain: acceptance rate and state of the chain. Default value to \code{n.iterations/100}. Set to \code{NULL} to avoid any info.
#' @param ... other arguments to be passed to function \code{log.likelihood}.
#' @inheritParams marginalLogLikelihoodDeterministic
#' @note By specifying \code{log.likelihood}, this function offers the possibility to run a MCMC independently of the method used for likelihood computation.
#' Although for a deterministic model \code{log.likelihood} is quite straightforward, fitting a stochastic model will require a call to a SMC algorithm to evaluate the log-likelihood.
#' Finally, this function can also implement an ABC procedure to decide whether to accept or reject the parameter set.
#'
#' The size of the proposal covariance matrix is adapted using the following formulae: \deqn{\Sigma_{n+1}=\sigma_n * \Sigma_n} with \eqn{\sigma_n=\sigma_{n-1}*exp(\alpha^n*(acc - 0.234))},
#' where \eqn{\alpha} is equal to \code{adapt.size.cooling} and \eqn{acc} is the acceptance rate of the chain.
#'
#' The shape of the proposal covariance matrix is adapted using the following formulae: \deqn{\Sigma_{n+1}=2.38^2/d * \Sigma_n} with \eqn{\Sigma_n} the empirical covariance matrix
#' and \eqn{d} is the number of estimated parameters in the model.
#' @references Roberts GO, Rosenthal JS. Examples of adaptive MCMC. Journal of Computational and Graphical Statistics. Taylor & Francis; 2009;18(2):349-67.
#' @export
#' @import tmvtnorm 
#' @importFrom lubridate as.period
#' @return a list with 3 elements:
#'\itemize{
#'	\item \code{trace} a \code{data.frame}. Each row contains a state of the chain (as returned by \code{target}) and the last column "weight" is the number of iteration the chain stayed in that state. This offers a compact storage for the full trace of the chain.
#'	\item \code{acceptance.rate} acceptance rate of the MCMC chain.
#'	\item \code{covmat.empirical} empirical covariance matrix of the target sample.
#'}
mcmcMH <- function(target, target.args, theta.init, gaussian.proposal=list(covmat=NULL, lower=NULL, upper=NULL), n.iterations, adapt.size.start=n.iterations, adapt.size.cooling=0.99, adapt.shape.start=n.iterations, print.info.every=n.iterations/100) {

	# initialise theta
	theta.current <- theta.init
	theta.propose <- theta.init

	# extract theta of gaussian proposal
	covmat.proposal <- gaussian.proposal$covmat
	lower.proposal <- gaussian.proposal$lower
	upper.proposal <- gaussian.proposal$upper

	# reorder all vector and matrix by names
	theta.names <- names(theta.init)
	covmat.proposal <- covmat.proposal[theta.names,theta.names]
	lower.proposal <- lower.proposal[theta.names]
	upper.proposal <- upper.proposal[theta.names]

	# covmat init
	covmat.proposal.init <- covmat.proposal
	start.adapt.size <- TRUE
	start.adapt.shape <- TRUE
	
	# find estimated theta
	theta.estimated.names <- names(which(diag(covmat.proposal)>0))

	# evaluate target at theta init
	target.theta.current <- do.call(target, c(list(theta=theta.current), target.args))

	# initialise trace data.frame
	trace <- data.frame(t(target.theta.current$trace), weight=1)

	# acceptance rate
	acceptance.rate <- 0

	# scaling factor for covmat size
	scaling.sd  <- 1

	# empirical covariance matrix (0 everywhere initially)
	covmat.empirical <- covmat.proposal
	covmat.empirical[,] <- 0

	# empirical mean vector
	theta.mean <- theta.current

	# if print.info.every is null never print info
	if(is.null(print.info.every)){
		print.info.every <- n.iterations + 1
	}

	start_iteration_time <- Sys.time()

	for (i.iteration in seq_len(n.iterations)) {

		# adaptive step
		if(i.iteration >= adapt.size.start && acceptance.rate*i.iteration < adapt.shape.start){
			if(start.adapt.size){
				message("\n---> Start adapting size of covariance matrix")
				start.adapt.size <- 0				
			}
			# adapt size of covmat until we get enough accepted jumps
			scaling.sd <- scaling.sd*exp(adapt.size.cooling^(i.iteration-adapt.size.start)*(acceptance.rate - 0.234))
			covmat.proposal <- scaling.sd^2*covmat.proposal.init

		}else if(acceptance.rate*i.iteration >= adapt.shape.start){
			if(start.adapt.shape){
				message("\n---> Start adapting shape of covariance matrix")
				# flush.console()
				start.adapt.shape <- 0
			}
			# adapt shape of covmat using optimal scaling factor for multivariate traget distributions
			covmat.proposal <- 2.38^2/length(theta.estimated.names)*covmat.empirical
		}

		# print info
		if(i.iteration%%round(print.info.every)==0){
			end_iteration_time <- Sys.time()
			x <- trace[nrow(trace),]
			x <- paste(paste0(names(x),"=",sprintf("%.2f",x)),collapse="|")
			suppressMessages(time.estimation <- round(as.period((end_iteration_time-start_iteration_time)*10000/round(print.info.every))))
			message("Iteration: ",i.iteration,"/",n.iterations," Time 10000 iter: ",time.estimation," Acceptance rate: ",sprintf("%.3f",acceptance.rate)," Scaling.sd: ",sprintf("%.3f",scaling.sd)," State:",x)
			start_iteration_time <- end_iteration_time
		}

		# propose another parameter set
		if(any(diag(covmat.proposal)[theta.estimated.names]<.Machine$double.eps)){
			print(covmat.proposal[theta.estimated.names,theta.estimated.names])
			stop("non-positive definite covmat",call.=FALSE)
		}
		theta.propose[theta.estimated.names] <- as.vector(rtmvnorm(1,mean=theta.current[theta.estimated.names],sigma=covmat.proposal[theta.estimated.names,theta.estimated.names],lower=lower.proposal[theta.estimated.names],upper=upper.proposal[theta.estimated.names]))

		# evaluate posterior of proposed parameter
		target.theta.propose <- do.call(target, c(list(theta=theta.propose), target.args))

		if(!is.finite(target.theta.propose$log.dist)){
			# if posterior is 0 then do not compute anything else and don't accept
			log.acceptance <- -Inf

		}else{

			# compute Metropolis-Hastings ratio (acceptance probability)
			log.acceptance <- target.theta.propose$log.dist - target.theta.current$log.dist + dtmvnorm(x=theta.current[theta.estimated.names],mean=theta.propose[theta.estimated.names],sigma=covmat.proposal[theta.estimated.names,theta.estimated.names],lower=lower.proposal[theta.estimated.names],upper=upper.proposal[theta.estimated.names],log=TRUE) - dtmvnorm(x=theta.propose[theta.estimated.names],mean=theta.current[theta.estimated.names],sigma=covmat.proposal[theta.estimated.names,theta.estimated.names],lower=lower.proposal[theta.estimated.names],upper=upper.proposal[theta.estimated.names],log=TRUE)

		}

		if(is.accepted <- (log(runif(1)) < log.acceptance)){
			# accept proposed parameter set
			trace <- rbind(trace,c(target.theta.propose$trace, weight=1))
			theta.current <- theta.propose
			target.theta.current <- target.theta.propose
		}else{
			# reject
			trace$weight[nrow(trace)] <- trace$weight[nrow(trace)] + 1
		}

		# update acceptance rate
		acceptance.rate <- acceptance.rate + (is.accepted - acceptance.rate)/i.iteration

		# update empirical covariance matrix
		tmp <- updateCovmat(covmat.empirical,theta.mean,theta.current,i.iteration)
		covmat.empirical <- tmp$covmat
		theta.mean <- tmp$theta.mean

	}

	return(list(trace=trace,acceptance.rate=acceptance.rate,covmat.empirical=covmat.empirical))
}



