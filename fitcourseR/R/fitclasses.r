
#'Constructor of fitparam object
#'
#'A \code{fitparam} object contains all the information necessary to fit a model parameter.
#' @param name character, name of the parameter.
#' @param value numeric, value of the parameter.
#' @param support numeric vector, range of the support of the parameter.
#' @param sd.proposal numeric, standard deviation of the proposal gaussian kernel.
#' @param prior list of 2 elements: 
#' \itemize{
#'	\item \code{distribution} character, name of density function for the prior distribution. IMPORTANT: it must be the name of recognized \R function (ex: "dnorm", "dunif").
#' 	\item \code{parameters} named vector of parameters for the prior distribution, as specified in the documentation of the \R function. For instance, the \code{mean} and \code{sd} of a normal distribution, or the \code{min} and \code{max} of a uniform distribution.
#' }
#' @export
#' @seealso getParameterValues, setParameterValues
#' @return a \code{fitparam} object
#' @examples \dontrun{
#' myParam <- fitparam(name="R0",value=3,support=c(0,Inf),sd.proposal=1,prior=list(distribution="uniform",parameters=c(min=1,max=100))) 
#'}
fitparam <- function(name=NULL,value=NULL,support=c(-Inf,Inf),sd.proposal=0,prior=list(distribution=NULL,parameters=NULL)) {

	if(is.null(name) || !is.character(name)){
		stop(sQuote("name")," argument must be a character",call.=FALSE)
	}
	if(is.null(value) || !is.numeric(value)){
		stop(sQuote("value")," argument must be numeric",call.=FALSE)
	}
	if(length(support)!=2 || !is.numeric(support)){
		stop(sQuote("support")," argument must be a numeric vector of length 2",call.=FALSE)
	}
	
	# reorder support just in case
	support <- sort(support)
	# test that default value is within the support range
	if(value < min(support) || value > max(support)){
		stop(sQuote("value")," argument is not within the support",call.=FALSE)
	}

	if(!is.numeric(sd.proposal) || sd.proposal<0 ){
		stop(sQuote("sd.proposal")," argument must be numeric and positive",call.=FALSE)
	}

	if(sd.proposal>0 && is.null(prior$distribution)){
		stop("Prior distribution must be provided (parameter is estimated)")
	}

	if(length(x <- setdiff(c("distribution","parameters"),names(prior)))){
		stop("element(s) ",paste(sQuote(x),collapse=", ")," is (are) missing in the prior list",call.=FALSE)
	}

	if(!is.null(prior$distribution)){
		# test that the prior is well specified
		value.prior <- do.call(prior$distribution,as.list(c(x=value,prior$parameters)))
		# test that value.prior is non-zero as it is the initial value for parameter exploration.
		if(value.prior==0){
			stop("The prior density is 0 for the value provided.",call.=FALSE)
		}
	}
	
	return(structure(list(
		name=name,
		value=value,
		support=support,
		sd.proposal=sd.proposal,
		prior=prior
		),class="fitparam"))
}

#'Check a list of fitparam objects
#'
#' @param list.fitparam named list of \code{\link{fitparam}} objects.
#' @export
checkListFitparam <- function(list.fitparam) {
	if(!all(x <- sapply(list.fitparam,inherits,what="fitparam"))){
		stop("elements ",paste(which(!x),collapse=",")," in list.fitparam argument are not fitparam objects")
	}
}

#'Get parameter values from a list of fitparam objects
#'
#'If \code{list.fitparam} is a list of \code{\link{fitparam}} objects, the value of the theta will be returned as a named vector, possibly restricting to estimated theta only.
#' @param list.fitparam named list of \code{\link{fitparam}} objects.
#' @export
#' @return a named vector of parameter values
#' @examples \dontrun{
#' list.fitparam <- list(fitparam(name="A",value=1,support=c(0,10),prior=list(distribution="uniform",theta=c(min=0,max=10))),fitparam(name="B",value=2),fitparam(name="C",value=2))
#' getParameterValues(list.fitparam)
#'}
getParameterValues <- function(list.fitparam) {

	checkListFitparam(list.fitparam)

	param.values <- sapply(list.fitparam,function(x) {x$value})
	names(param.values) <- sapply(list.fitparam,function(x) {x$name})

	return(param.values)
}


#'Replace specified parameter values with new values, in a list of fitparam objects.
#'
#'If \code{list.fitparam} is a list of \code{\link{fitparam}} objects, the value of the theta will be replaced by those of \code{new.values}, by matching parameter names with names.
#' @param list.fitparam list of \code{\link{fitparam}} objects.
#' @param new.values named vector containing the new parameter values; names must match those of \code{list.fitparam}.
#' @export
#' @seealso fitparam
#' @return a named list of \code{\link{fitparam}} objects.
#' @examples \dontrun{
#' list.fitparam <- list(fitparam(name="A",value=1),fitparam(name="B",value=2),fitparam(name="C",value=2))
#' list.fitparam <- setParameterValues(list.fitparam,new.value=c("A"=10,"B"=5))
#' print(list.fitparam)
#'}
setParameterValues <- function(list.fitparam,new.values) {

	checkListFitparam(list.fitparam)
	
	# name the list
	names(list.fitparam) <- sapply(list.fitparam,function(x) {x$name})

	for(x in names(new.values)){
		list.fitparam[[x]]$value <- new.values[[x]]
	}

	# uname the list
	names(list.fitparam) <- NULL

	return(list.fitparam)

}

#'Constructor of fitmodel object
#'
#'A \code{fitmodel} object contains all the information necessary to simulate and fit a model during the course. When a model is created, the constructor performs a serie of checks
#'on the arguments provided by the user in order to make sure that the different elements of the \code{fitmodel} will be compatible both with the functions coded during the course and the functions 
#'available in the \code{fitcourseR} package. The latter can be used as a correction. 
#' @param name character, name of the model (required)
#' @param state.variables character vector, names of the state variables i.e. \code{c("S","I","R","Incidence")} (required)
#' @param list.fitparam list of \code{\link{fitparam}} objects. This list will be used to create a named vector of model parameters (\code{theta}), a function to compute the prior (\code{log.prior}) and a list of parameters (including a covariance matrix) for the gaussian proposal kernel (\code{gaussian.proposal}). (required)
#' @param initialise.state \R function to initialise the state of the models using the model parameters \code{theta} (optional). This function takes one argument:
#' \itemize{
#' 	\item \code{theta} named numeric vector of model parameters 
#' } 
#' and returns a numeric vector of initial states, whose names match the vector \code{state.variables}.
#' @param simulate.model \R function to simulate forward the model. This function takes 3 arguments:
#' \itemize{
#' \item \code{theta} named vector of model parameters.
#' \item \code{state.init} named vector of initial state of the model.
#' \item \code{times} time sequence for which state of the model is wanted; the first value of times must be the initial time.
#' } 
#' and returns a \code{data.fame} containing the values of the state variables (1 per column) at each observation time (1 per row). 
#' @param generate.observation \R-function to generate simulated data from a simulated trajectory using an observation model. This function takes 2 arguments:
#' \itemize{
#' \item \code{model.traj} data.frame of simulated trajectories, as returned by \code{simulate.model}.
#' \item \code{theta} named vector of model parameters.
#' }
#' and return the \code{model.traj} data.frame with an additional "observation" column generated by the observation model. (optional)
#' @param log.prior.fitparam \R function to evaluate the log-prior distribution, using the information about priors contained in \code{list.fitparam} (see note below). This function takes one argument
#' \itemize{
#' 	\item \code{list.fitparam} list of \code{\link{fitparam}} objects.
#' }
#' and returns the logged value of the prior density distribution. Required if at least one parameter is estimated. See note below.
#' @param data \code{data.frame} of times and observations the model should be fitted to. \code{data} must contain a column named "time" with all observation times. (optional)
#' @param log.likelihood \R function to evaluate the log-likelihood of the data given a model parameter \code{theta} and a simulated trajectory \code{model.traj}. This function takes 3 arguments:
#' \itemize{
#' \item \code{data} data.frame containing all or part of the observations.
#' \item \code{theta} named vector of model parameters.
#' \item \code{model.traj} data.frame containing the state of the model at the successive observation times, as returned by \code{simulate.model}.
#' }
#' and return the log-likelihood. (optional)
#' @param distance.ABC \R function to evaluate one or more summary distances between the data and a simulated observation. This function takes 2 arguments:
#' \itemize{
#' \item \code{data} data.frame containing all or part of the observations.
#' \item \code{model.traj.obs} data.frame containing simulated observations, as returned by \code{generate.observation}.
#' }
#' and return a numeric vector of summary distances. (optional)
#' @param gaussian.proposal list of parameters for the - potentially truncated - gaussian proposal distribution of the MCMC. Contains 3 elements:
#'\itemize{
#'	\item \code{covmat} covariance matrix. Must have named rows and columns with at least all estimated theta. (optional)
#'	\item \code{lower} named vector of lower truncation points. Will be set to \code{-Inf} by default.
#'	\item \code{upper} named vector of upper truncation points. Will be set to \code{Inf} by default.
#'} 
#' @param verbose if \code{TRUE}, print details of the test performed to check validity of the arguments
#' @export
#' @seealso fitparam compositeLogPrior
#' @note \itemize{
#' \item The \code{log.prior.fitparam} function can take advantage of the way univariate prior distributions are defined in the \code{\link{fitparam}} object to compute a composite log-prior. See \code{\link{compositeLogPrior}} for more details.
#' }
#' @return a \code{fitmodel} object that is a \code{list} of 9 elements:
#' \itemize{
#' 	\item \code{state.variables} vector names of the state variables.
#' 	\item \code{theta} named vector of the values of model parameters.
#' 	\item \code{initialise.state} \R function to initialise the states of the model; usage: \code{initialise.state(theta)}.
#' 	\item \code{simulate.model} \R function to simulate forward the model; usage: \code{simulate.model(theta,state.init,times)}.
#' 	\item \code{generate.observation} \R function to generate simulated observation; usage: \code{generate.observation(model.traj, theta)}.
#' 	\item \code{log.prior} \R function to evaluate the log-prior; usage: \code{log.prior(theta)}.
#' 	\item \code{data} data.frame with observation times and measures.
#' 	\item \code{log.likelihood} \R function to evaluate the log-likelihood; usage: \code{log.likelihood(data, theta, model.traj)}.
#' 	\item \code{gaussian.proposal} parameters of the gaussian proposal kernel. A list of 3 elements:
#' 	\itemize{
#' 		\item \code{covmat} covariance matrix
#' 		\item \code{lower} lower truncation vector
#' 		\item \code{upper} upper truncation vector
#' 	}
#' } 
#' @example inst/examples/example-fitmodel.r
fitmodel <- function(name=NULL,state.variables=NULL, list.fitparam=NULL, initialise.state=NULL, simulate.model=NULL, generate.observation=NULL, log.prior.fitparam=NULL, data=NULL, log.likelihood=NULL, distance.ABC=NULL ,gaussian.proposal=list(covmat=NULL,lower=NULL,upper=NULL), verbose=TRUE){
	
	# mandatory
	if(!is.character(name)){
		stop(sQuote("name")," argument is not a character")
	}
	if(!is.character(state.variables)){
		stop(sQuote("state.variables")," argument is not a character vector")
	}
	
	checkListFitparam(list.fitparam)

	if(!is.null(simulate.model)) {
		if(!is.function(simulate.model)){
			stop(sQuote("simulate.model")," argument is not an R function")
		}
		if(is.null(initialise.state)) {
			stop(sQuote("initialise.state")," must be provided when", sQuote("simulate.model"), " is given")
		} else if(!is.function(initialise.state)) {
			stop(sQuote("initialise.state")," argument is not an R function")
		}
	}
	# optional (required only if at least one parameter is estimated)
	if(!is.null(log.prior.fitparam) && !is.function(log.prior.fitparam)){
		stop(sQuote("log.prior.fitparam")," argument is not an R function")
	}
	if(!is.null(data) && !is.data.frame(data)){
		stop(sQuote("data")," argument is not a data.frame")
	}
	if(!is.null(log.likelihood) && !is.function(log.likelihood)){
		stop(sQuote("log.likelihood") ," argument is not an R function")
	}
	if(!is.null(distance.ABC) && !is.function(distance.ABC)){
		stop(sQuote("distance.ABC")," argument is not an R function")
	}
	if(!is.null(gaussian.proposal$covmat) && !is.matrix(gaussian.proposal$covmat)){
		stop(sQuote("covmat"), " element in ",sQuote("gaussian.proposal") ," argument is not a matrix")
	}
	if(!is.null(gaussian.proposal$lower) && !is.vector(gaussian.proposal$lower)){
		stop(sQuote("lower"), " element in ",sQuote("gaussian.proposal") ," argument is not a vector")
	}
	if(!is.null(gaussian.proposal$upper) && !is.vector(gaussian.proposal$upper)){
		stop(sQuote("upper"), " element in ",sQuote("gaussian.proposal") ," argument is not a vector")
	}
	if(!is.null(generate.observation) && !is.function(generate.observation)){
		stop(sQuote("generate.observation")," argument is not an R function")
	}

	# create the named vector of theta
	theta <- getParameterValues(list.fitparam)
	names.theta <- names(theta)
	# name list.fitparam for convenience	
	names(list.fitparam) <- names.theta

	# check initialise.state
	if(!is.null(initialise.state)) {
		if(verbose){
			cat("--- check initialise.state\n")
		}
		fun_args <- c("theta")
		if(!all(x <- fun_args%in%names(formals(initialise.state)))){
			stop("argument(s) ",paste(fun_args[!x],collapse=", ")," missing in function initialise.state, see documentation.")
		}

		test.initial.state <- initialise.state(theta=theta)

		if(verbose){
			cat("initialise.state(theta) should return a non-negative numeric vector of dimension ",length(state.variables)," with names: ",paste(state.variables,collapse=", "),"\nTest:\n")
			print(test.initial.state)
		}
		if(!is.vector(test.initial.state)){
			stop("initialise.state must return a vector")
		}
		if(!all(x <- state.variables%in%names(test.initial.state))){
			stop("State variable(s) missing in the vector returned by the function initialise.state:",paste(state.variables[!x],collapse=", "))
		}
		if(!all(x <- names(test.initial.state)%in%state.variables)){
			stop("The following state variables should not be returned by the function initialise.state since they are not defined in state.variables:",paste(names(test.initial.state)[!x],collapse=", "))
		}
		if(any(x <- test.initial.state<0)){
			stop("The following state variables are negative:",paste(test.initial.state[x],collapse=", "))
		}
		if(verbose){
			cat("--> initialise.state looks good!\n")
		}
	}

	# check simulate.model
	if(!is.null(simulate.model)) {
		if(verbose){
			cat("--- check simulate.model\n")
		}
                                        # check arguments
		fun_args <- c("theta","state.init","times")
		if(!(all(x <- fun_args%in%names(formals(simulate.model))))){
			stop("argument(s) ",paste(fun_args[!x],collapse=", ")," missing in function simulate.model, see documentation.")
		}

		times <- 0:10
		test.traj <- simulate.model(theta=theta,state.init=test.initial.state,times=times)
                                        # must return a data.frame of dimension 11x(length(state.variables)+1)
		if(verbose){
			cat("simulate.model(theta, test.initial.state, times=0:10) should return a non-negative data.frame of dimension",length(times),"x",length(state.variables)+1,"with column names:",paste(c("time",state.variables),collapse=", "),"\nTest:\n")
			print(test.traj)
		}
		if(!is.data.frame(test.traj)){
			stop("simulate.model must return a data.frame")
		}
		if(!all(x <- c("time",state.variables)%in%names(test.traj))){
			stop("Column(s) missing in the data.frame returned by simulate.model:",paste(c("time",state.variables)[x],collapse=", "))
		}
		if(!all(x <- names(test.traj)%in%c("time",state.variables))){
			warning("The following columns are not required in the data.frame returned by simulate.model:",paste(names(test.traj)[!x],collapse=", "))
		}
		if(any(test.traj$time!=times)){
			stop("The time column of the data.frame returned by simulate.model is different from its times argument",call.=FALSE)
		}
		if(any(test.traj<0)){
			stop("simulate.model returned negative values during the test, use verbose argument of fitmodel to check")
		}
		if(verbose){
			cat("--> simulate.model looks good!\n")
		}
	}

	# check generate.observation
	if(!is.null(generate.observation)) {
	    # check arguments
		fun_args <- c("model.traj","theta")
		if(!(all(x <- fun_args%in%names(formals(generate.observation))))){
			stop("argument(s) ",paste(fun_args[!x],collapse=", ")," missing in function generate.observation, see documentation.")
		}

		test.generate.observation <- generate.observation(test.traj, theta)
		if(verbose){
			cat("generate.observation(test.traj, theta) should return a non-negative data.frame of dimension",nrow(test.traj),"x",ncol(test.traj)+1,"with column names:",paste(c(names(test.traj),"observation"),collapse=", "),"\nTest:\n")
			print(test.generate.observation)
		}
		if(!is.data.frame(test.generate.observation)){
			stop("generate.observation must return a data.frame")
		}
		if(!all(x <- c(names(test.generate.observation),"observation")%in%names(test.generate.observation))){
			stop("Column(s) missing in the data.frame returned by generate.observation:",paste(c("time",state.variables)[x],collapse=", "))
		}
		if(!all(x <- names(test.generate.observation)%in%c(names(test.generate.observation),"observation"))){
			warning("The following columns are not required in the data.frame returned by generate.observation:",paste(names(test.generate.observation)[x],collapse=", "))
		}
		if(nrow(test.generate.observation)!=nrow(test.traj)){
			stop("The data.frame returned by generate.observation must have the same number of rows as the model.traj argument",call.=FALSE)
		}
		if(any(test.generate.observation$observation<0)){
			stop("generate.observation returned negative observation during the test, use verbose argument of fitmodel to check")
		}
		if(verbose){
			cat("--> generate.observation looks good!\n")
		}
	}

	# build covariance matrix
	covmat.proposal <- gaussian.proposal$covmat

	if(!is.null(covmat.proposal)){

		# check covmat proposal
		if(verbose){
			cat("Check covariance matrix of the gaussian proposal kernel:\n")
			print(covmat.proposal)
		}
		# must be a square matrix with same row and col names
		if(nrow(covmat.proposal)!=ncol(covmat.proposal)){
			stop("gaussian.proposal$covmat is not a square matrix")
		}
		if(length(setdiff(rownames(covmat.proposal),colnames(covmat.proposal)))){
			stop("row names and col names of gaussian.proposal$covmat don't match")
		}

		# check whether all names correspond to existing theta
		names.unknown.theta <- setdiff(rownames(covmat.proposal),names.theta)
		if(length(names.unknown.theta)){
			stop("The following parameter names in gaussian.proposal$covmat do not match those of list.fitparam: ",paste(names.unknown.theta,collapse=", "))
		}

		# check that all theta are in gaussian.proposal$covmat
		names.missing.theta <- setdiff(names.theta,rownames(covmat.proposal))
		if(length(names.missing.theta)){
			# if missing theta => set their variances to 0 (not estimated)
			warning("Parameter(s): ",paste(names.missing.theta,collapse=", ")," missing in gaussian.proposal$covmat. The matrix will be updated to include missing theta but they won't be estimated (proposal variance = 0).",call.=FALSE)
			new.covmat.proposal <- matrix(0,ncol=length(names.theta),nrow=length(names.theta),dimnames=list(names.theta,names.theta))
			new.covmat.proposal[rownames(covmat.proposal),colnames(covmat.proposal)] <- covmat.proposal
			covmat.proposal <- new.covmat.proposal
		}
		if(verbose){
			cat("--> gaussian.proposal$covmat looks good!\n")
		}


	}else{
		# use list.fitparam
		covmat.proposal <- diag(sapply(list.fitparam,function(x) {(x$sd.proposal)^2}), nrow = length(list.fitparam))
		rownames(covmat.proposal) <- colnames(covmat.proposal) <- names.theta

	}	

	gaussian.proposal$covmat <- covmat.proposal

	# extract list of estimated theta
	estimated.theta <- names(which(diag(covmat.proposal)>0))

	# if any estimated theta do some tests on required arguments
	if(length(estimated.theta)){

		# check that all estimated theta have a prior (if so, prior has alreday been tested in fitparam)
		if(any(x <- sapply(list.fitparam[estimated.theta],function(x) {is.null(x$prior$distribution)}))){
			stop("Prior argument in fitparam must be defined for the following theta:",paste(estimated.theta[x],collapse=", "),call.=FALSE)
		}

		if(is.null(log.prior.fitparam)){
			stop(sQuote("log.prior.fitparam")," argument must be provided because at least 1 parameter is marked as estimated")
		}
		if(is.null(data)){
			stop(sQuote("data")," argument must be provided because at least 1 parameter is marked as estimated")
		}
		if(is.null(log.likelihood) && is.null(distance.ABC)){
			stop("Either",sQuote("log.likelihood")," or ", sQuote("distance.ABC")," arguments must be provided because at least 1 parameter is estimated")
		}

		if(!is.null(lower <- gaussian.proposal$lower)){
			## check lower truncations

			# subset lower to parameter names
			lower <- lower[intersect(names(lower),names.theta)]
			# fill missing theta with -Inf
			names.default <- setdiff(names.theta,names(lower))
			lower.default <- rep(-Inf,length(names.default))
			names(lower.default) <- names.default
			# bind
			lower <- c(lower,lower.default)
			# reorder
			lower <- lower[names.theta]
		}else{
			#use list.fitparam
			lower <- sapply(list.fitparam,function(x) {min(x$support)})
			names(lower) <- names.theta
		}
		
		if(!is.null(upper <- gaussian.proposal$upper)){
			## check upper truncations
			upper <- upper[intersect(names(upper),names.theta)]
			# fill missing theta with -Inf
			names.default <- setdiff(names.theta,names(upper))
			upper.default <- rep(Inf,length(names.default))
			names(upper.default) <- names.default
			# bind
			upper <- c(upper,upper.default)
			# reorder
			upper <- upper[names.theta]
		}else{
			#use list.fitparam
			upper <- sapply(list.fitparam,function(x) {max(x$support)})
			names(upper) <- names.theta
		}		

		# check lower < upper
		if(any(x <- (lower >= upper))){
			stop("lower must be < than upper for parameter(s): ",paste(names(x)[x],collapse=", "),call.=FALSE)
		}
		gaussian.proposal$lower <- lower
		gaussian.proposal$upper <- upper

		if(verbose){
			cat("--> Gaussian proposal will be truncated to:\n")
			print(data.frame(lower=lower,upper=upper))
		}

	}

	if (is.null(log.prior.fitparam)) {
		log.prior <- NULL
	} else {

		# check arguments
		fun_args <- c("list.fitparam")
		if(!(all(x <- fun_args%in%names(formals(log.prior.fitparam))))){
			stop("arguments ",paste(fun_args[!x],collapse=", ")," missing in function log.prior.fitparam, see documentation.")
		}

		log.prior <- function(theta) {
			return(log.prior.fitparam(setParameterValues(list.fitparam,theta)))
		}

		# test it
		test.log.prior <- log.prior(theta)
		if(verbose){
			cat("log.prior(theta) should return a single finite value\nTest:",test.log.prior,"\n")
		}
		if(!(!is.na(test.log.prior) && (is.finite(test.log.prior)))){
			stop("log.prior must return a finite value for initial parameter values")
		}
		if(verbose){
			cat("--> log.prior looks good!\n")
		}
	}

	# data must have a column named time, should not start at 0
	if (!is.null(data)) {
		if(!"time"%in%names(data)){
			stop("data argument must have a column named \"time\"")
		}else if(data$time[1]==0){
			stop("the first observation time in data argument should not be 0")
		}
	}

    # check log.likelihood
	# check arguments, return value
	if (!is.null(log.likelihood)) {

		if(is.null(simulate.model)) {
			stop(sQuote("simulate.model")," must be provided when", sQuote("log.likelihood"), " is given")
		} else if(!is.function(simulate.model)) {
			stop(sQuote("simulate.model")," argument is not an R function")
		}

		if (is.null(data)) {
			stop("data argument must be provided for log.likelihood function")
		}
		# check arguments
		fun_args <- c("data","theta","model.traj")
		if(!(all(x <- fun_args%in%names(formals(log.likelihood))))){
			stop("argument(s) ",paste(fun_args[!x],collapse=", ")," missing in function log.likelihood, see documentation.")
		}

        # test it
		test.log.likelihood <- log.likelihood(data=data,theta=theta,model.traj=test.traj)

		if(verbose){
			cat("log.likelihood(model.traj,data,theta) should return a single value\nTest:",test.log.likelihood,"\n")
		}
		if(is.na(test.log.likelihood) || (test.log.likelihood > 0)){
			stop("log.likelihood must return a non-positive value")
		}
		if(verbose){
			cat("--> log.likelihood looks good!\n")
		}
	}


	 # check distance.ABC
	# check arguments, return value
	if (!is.null(distance.ABC)) {

		if (is.null(data)) {
			stop(sQuote("data") ," argument must be provided for ABC inference")
		}
		if (is.null(generate.observation)){
			stop(sQuote("generate.observation")," argument must be provided for ABC inference")
		}
		# check arguments
		fun_args <- c("model.traj.obs","data")
		if(!(all(x <- fun_args%in%names(formals(distance.ABC))))){
			stop("argument(s) ",paste(fun_args[!x],collapse=", ")," missing in function ",sQuote("distance.ABC")," see documentation.")
		}

        # test it
        # test.generate.observation has been succesfully generated above (otherwise generate.observation is missing or corrupt and fitmodel has already thrown an error)
		test.distance.ABC <- distance.ABC(data=data,model.traj.obs=test.generate.observation)

		if(verbose){
			cat("distance.ABC(model.traj.obs,data) should return a numeric vector\nTest:",test.distance.ABC,"\n")
		}
		if(is.na(test.distance.ABC) || !is.numeric(test.distance.ABC)){
			stop("distance.ABC must return a numerci vector")
		}
		if(verbose){
			cat("--> distance.ABC looks good!\n")
		}
	}

	# create and return object
	return(structure(list(
		name=name, 
		state.variables=state.variables,
		theta=theta,
		initialise.state=initialise.state,
		simulate.model=simulate.model,
		generate.observation=generate.observation,
		log.prior=log.prior,
		data=data,
		log.likelihood=log.likelihood,
		distance.ABC=distance.ABC,
		gaussian.proposal=gaussian.proposal
		),class="fitmodel"))

}
