#' Initialise SEIT2L model
#'
#' Set the initial value of the state space of the SEIT2L model for a given set of theta. 
#' @param theta named vector of parameters for the SEIT2L model.
#' @note This function throws an error if the sum of all compartments is higher than the population size N.
#' @export
#' @seealso fitmodel, generateSEIT2LmodelTDC
#' @return a named vector of initial conditions for the SEIT2L model.
SEIT2L_initialiseState <- function(theta) {

	# constant pop size
	N <- theta[["N"]]

	# number of infected and immune
	I <- round(theta[["pI0"]]*N)
	L <- round(theta[["pL0"]]*N)

	if(I+L>N){
		stop("Initial conditions not valid")
	}

	return(c(S=N-I-L,E=0,I=I,T1=0,T2=0,L=L,Inc=0))
}


#' Log-prior for SEIT2L model
#'
#' Compute the log of the prior density distribution for the parameter set \code{theta}, by composing the univariate priors of all estimated theta.
#' @param list.fitparam list of \code{fitparam} objects.
#' @note This function returns \code{-Inf} if the sum of all compartments is higher than the population size N.
#' @export
#' @seealso fitmodel, compositeLogPrior
#' @return the value of the log-prior.
SEIT2L_logPrior <- function(list.fitparam) {

	theta <- getParameterValues(list.fitparam)

	# check constraint, due to the constant population size, to avoid negative proportion of susceptible as initial condition
	if((theta[["pI0"]]+theta[["pL0"]]) > 1){
		return(-Inf)
	}

	return(compositeLogPrior(list.fitparam))

}


#'Deterministic simulation of SEIT2L model
#'
#'Solves the system of ordinary differential equations for the SEIT2L model using the \code{\link[deSolve]{ode}} function.
#' @param state.init a named vector with the initial state values for the SEIT2L model. 
#' @param times time sequence for which output is wanted; the first value of \code{times} must be the initial time.
#' @inheritParams SEIT2L_initialiseState
#' @export
#' @import deSolve
#' @return a \code{data.fame} containing the values of the state variables (columns) at each observation times (rows). 
SEIT2L_simulateDeterministic <- function(theta,state.init,times) {

	SEIT2L_ode <- function(time, state, theta) {

		# param
		beta <- theta[["R0"]]/theta[["IP"]]
		N <- theta[["N"]]
		epsilon <- 1/theta[["LP"]]
		nu <- 1/theta[["IP"]]
		alpha <- theta[["alpha"]]
		tau <- 2/theta[["TIP"]]

		# states
		S <- state[["S"]]
		E <- state[["E"]]
		I <- state[["I"]]
		T1 <- state[["T1"]]
		T2 <- state[["T2"]]
		L <- state[["L"]]
		Inc <- state[["Inc"]]

		dS <- -beta*S*I/N + (1-alpha)*tau*T2
		dE <- beta*S*I/N - epsilon*E
		dI <- epsilon*E - nu*I
		dT1 <- nu*I - tau*T1
		dT2 <- tau*T1 - tau*T2
		dL <- alpha*tau*T2
		dInc <- epsilon*E

		return(list(c(dS,dE,dI,dT1,dT2,dL,dInc)))
	}


	traj <- as.data.frame(ode(state.init, times, SEIT2L_ode, theta))

	return(traj)

}

#'Stochastic simulation of SEIT2L model
#'
#'Simulate realisation of the stochastic version of the SEIT2L model using the \code{simulateModelStochastic} function.
#' @inheritParams SEIT2L_simulateDeterministic
#' @export
#' @seealso fitmodel, simulateModelStochastic
#' @return a \code{data.fame} containing the values of the state variables (columns) at each observation times (rows). 
SEIT2L_simulateStochastic <- function(theta,state.init,times) {

	SEIT2L_transitions <- list(
		c(S=-1,E=1),# infection
		c(E=-1,I=1,Inc=1),# infectiousness + incidence
		c(I=-1,T1=1),# recovery + short term protection
		c(T1=-1,T2=1),# recovery + short term protection
		c(T2=-1,L=1),# efficient long term protection
		c(T2=-1,S=1)# deficient long term protection
		)

	SEIT2L_rateFunc <- function(state,theta,t) {

		# param
		beta <- theta[["R0"]]/theta[["IP"]]
		N <- theta[["N"]]
		epsilon <- 1/theta[["LP"]]
		nu <- 1/theta[["IP"]]
		alpha <- theta[["alpha"]]
		tau <- 2/theta[["TIP"]]

		# states
		S <- state[["S"]]
		E <- state[["E"]]
		I <- state[["I"]]
		T1 <- state[["T1"]]
		T2 <- state[["T2"]]
		L <- state[["L"]]
		Inc <- state[["Inc"]]

		return(c(
			beta*S*I/N, # infection
			epsilon*E, # infectiousness + incidence
			nu*I, # recovery + short term protection
			tau*T1, # progression short term protection
			alpha*tau*T2, # efficient long term protection
			(1-alpha)*tau*T2 # deficient long term protection
			)
		)
	}

	return(simulateModelStochastic(theta,state.init,times,SEIT2L_transitions,SEIT2L_rateFunc))

}


#'Likelihood of the data for SEIT2L model
#'
#'Computes the log-likelihood of a subset of the data for a fixed trajectory and under a poisson observation process.
#' @param data subset of the \code{\link{FluTdC1971}} dataset that contains the observations corresponding to \code{model.traj}.
#' @param model.traj simulated trajectory, as returned by \code{\link{SEIT2L_simulateDeterministic}} or \code{\link{SEIT2L_simulateStochastic}}.
#' @inheritParams SEIT2L_initialiseState
#' @note This function can be used to compute the likelihood of the data given a deterministic trajectory or the weight of a SMC particle at a given observation time.
#' @export
#' @seealso SEIT2L_generateObservation
#' @return the log-likelihood value.
SEIT2L_logLikelihood <- function(data, theta, model.traj){

	# daily incidence needed
	daily.incidence <- diff(model.traj$Inc)

	# keep only data incidence corresponding to simulated times
	data <- subset(data,time%in%model.traj$time[-1]) # [-1] to remove initial simulation time

	x <- sum(dpois(x=data$Inc,lambda=theta[["rho"]]*daily.incidence,log=TRUE))

	return(x)
}


#'Generate an observed incidence time series
#'
#'Generate a daily incidence time serie under a Poisson observation process.  
#' @inheritParams SEIT2L_logLikelihood 
#' @inheritParams SEIT2L_initialiseState
#' @export
#' @seealso SEIT2L_likelihood
#' @return the \code{model.traj} data.frame with an additional variable: "observation".
SEIT2L_generateObservation <- function(model.traj, theta){

	# daily incidence needed
	daily.incidence <- diff(model.traj$Inc)

	x <- rpois(length(daily.incidence),lambda=theta[["rho"]]*daily.incidence)

	model.traj$observation <- c(0,x)

	return(model.traj)
}


#'ABC distance with oscillations
#'
#'This positive distance is the mean squared differences between the simulation and the observation, divided by the square of the number of times the simulation oscillates around the observation.
#' @param model.traj.obs \code{data.frame} of simulated trajectory with observation, as returned by \code{\link{SEIT2L_generateObservation}}.
#' @param data \code{data.frame} of times and observations. Must have two columns: \code{time} and \code{Inc}.
#' @export
#' @seealso distanceOscillation
#' @examples \dontrun{
#' # Suppose we observed a time series:
#' data <- data.frame(time=1:7,Inc=c(1,3,5,7,5,3,1))
#' # and we have two simulated time series:
#' traj1 <- data.frame(time=1:7,observation=c(3,5,7,9,7,5,3))
#' traj2 <- data.frame(time=1:7,observation=c(3,5,3,5,7,5,3))
#' # traj1 is consistently above data and traj2 oscillates around data:
#' plot(data$time,data$Inc,t='l',ylim=c(0,10))
#' lines(traj1$time,traj1$observation,col="red")
#' lines(traj2$time,traj2$observation,col="blue")
#' # While the squared differences are the same, we obtain a smaller distance for traj2:
#' d1 <- SEIT2L_distanceOscillation(traj1,data)
#' # d1 = 4
#' d2 <- SEIT2L_distanceOscillation(traj2,data)
#' # d2 = 1.3
#'}
SEIT2L_distanceOscillation <- function(model.traj.obs, data) {

	# match model and data on time
	keep.time <- intersect(model.traj.obs$time,data$time)
	model.traj.obs <- subset(model.traj.obs,time%in%keep.time)
	data <- subset(data,time%in%keep.time)
	
	x <- model.traj.obs$observation
	y <- data$Inc

	return(distanceOscillation(x,y))
}


#'Create the SEIT2L model as a fitmodel object
#'
#'This function returns a fitmodel object contaning the (either deterministic or stochastic) SEIT2L model, to be fitted to the \code{\link{FluTdC1971}} dataset.
#' @param deterministic if \code{TRUE} then the model will be deterministic, and stochastic otherwise.
#' @inheritParams fitmodel
#' @export
#' @import plyr
#' @return a fitmodel object
SEIT2L_createModelTdC <- function(deterministic=TRUE, verbose=TRUE) {

	# define theta using fitparam
	R0 <- fitparam(name="R0",value=10,support=c(0,Inf),sd.proposal=1,prior=list(distribution="dunif",parameters=c(min=1,max=100))) 

	LatentPeriod <- fitparam(name="LP",value=2,support=c(0,Inf),sd.proposal=0.5,prior=list(distribution="dunif",parameters=c(min=0,max=7))) 

	InfectiousPeriod <- fitparam(name="IP",value=3,support=c(0,Inf),sd.proposal=0.5,prior=list(distribution="dunif",parameters=c(min=0,max=30))) 

	TemporaryImmunePeriod  <- fitparam(name="TIP",value=10,support=c(0,Inf),sd.proposal=2,prior=list(distribution="dunif",parameters=c(min=0,max=50))) 

	ProbLongTermImmunity <- fitparam(name="alpha",value=0.5,support=c(0,1),sd.proposal=0.1,prior=list(distribution="dunif",parameters=c(min=0,max=1))) 

	ReportingRate <- fitparam(name="rho",value=0.7,support=c(0,2),sd.proposal=0.1,prior=list(distribution="dunif",parameters=c(min=0,max=2))) 

	proportionI0 <- fitparam(name="pI0",value=2/284,support=c(0,1),sd.proposal=1/284,prior=list(distribution="dunif",parameters=c(min=1/284,max=5/284))) 
	
	proportionL0 <- fitparam(name="pL0",value=0.1,support=c(0,1),sd.proposal=0.01,prior=list(distribution="dunif",parameters=c(min=0.0,max=0.5))) 
	
	PopSize <- fitparam(name="N",value=284) 

	# load and rename data
	data("FluTdC1971",envir = environment())
	data <- rename(FluTdC1971,c("day"="time","incidence"="Inc"))[c("time","Inc")]

	# simulator
	if(deterministic){
		simulate.model <- SEIT2L_simulateDeterministic
	}else{
		simulate.model <- SEIT2L_simulateStochastic
	}

	# create fitmodel
	SEIT2L <- fitmodel(
		verbose=verbose,
		name="SEIT2L",
		state.variables=c("S","E","I","T1","T2","L","Inc"),
		list.fitparam=list(R0,LatentPeriod,InfectiousPeriod,TemporaryImmunePeriod,ProbLongTermImmunity,ReportingRate,proportionI0,proportionL0,PopSize), 
		initialise.state=SEIT2L_initialiseState,
		log.prior.fitparam=SEIT2L_logPrior,
		simulate.model=simulate.model, 
		generate.observation=SEIT2L_generateObservation, 
		data=data, 
		log.likelihood=SEIT2L_logLikelihood,
		distance.ABC=SEIT2L_distanceOscillation
		) 

	return(SEIT2L)
}
