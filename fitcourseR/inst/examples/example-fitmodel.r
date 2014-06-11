# create a simple stochastic SIR model with constant population size

N <- 7.3e+6
# define model parameters using the fitparam class
R0 <- fitparam(name="R0",value=15,support=c(0,Inf),sd.proposal=1, prior=list(distribution="dunif",parameters=c(min=1,max=100)))

InfectiousPeriod <- fitparam(name="IP",value=2,support=c(0,Inf),sd.proposal=0.5, prior=list(distribution="dunif",parameters=c(min=0,max=30)))

ReportingRate <- fitparam(name="rho",value=0.7,support=c(0,2),sd.proposal=0.1, prior=list(distribution="dunif",parameters=c(min=0,max=2)))

proportionI0 <- fitparam(name="pI0",value=30/N,support=c(0,1),sd.proposal=1/N, prior=list(distribution="dunif",parameters=c(min=1/N,max=300/N)))
proportionR0 <- fitparam(name="pR0",value=0.9,support=c(0,1),sd.proposal=0.01, prior=list(distribution="dunif",parameters=c(min=0.5,max=1)))

PopSize <- fitparam(name="N",value=7.3e+6)

# function to initialise the model
SIR_initialiseState <- function(theta) {

	# constant pop size
	N <- theta[["N"]]

	# number of infected and immune
	I <- round(theta[["pI0"]]*N)
	R <- round(theta[["pR0"]]*N)

	if(I+R>N){
		stop("Initial conditions not valid")
	}

	return(c(S=N-I-R,I=I,R=R,Inc=0))
}



SIR_simulateDeterministic <- function(theta,state.init,times) {

	SIR_ode <- function(time,state,theta) {

        ## param
		beta <- theta[["R0"]]/theta[["IP"]]
		gamma <- 1/theta[["IP"]]
		N <- theta[["N"]]

        ## states
		S <- state["S"]
		I <- state["I"]
		R <- state["R"]
		Inc <- state["Inc"]

		dS <- -beta*S*I/N
		dI <- beta*S*I/N-gamma*I
		dR <- gamma*I
		Inc <- beta*S*I/N

		return(list(c(dS,dI,dR,Inc)))
	}

	trajectory <- data.frame(ode(y=state.init,times=times,func=SIR_ode,parms=theta))

	return(trajectory)
}

SIR_simulateStochastic <- function(theta,state.init,times) {

	# transitions
	SIR_transitions <- list(
		c(S=-1,I=1,Inc=1),# infectiousness + incidence
		c(I=-1,R=1)# recovery
		)

    # rates
	SIR_rateFunc <- function(x,theta,t) {

		# extract theta and states for convenience	when writing the rates
		beta <- theta[["R0"]]/theta[["IP"]]
		nu <- 1/theta[["IP"]]
		N <- theta[["N"]]
		S <- x["S"]
		I <- x["I"]

		return(c(
			beta*S*I/N, # infectiousness + incidence
			nu*I # recovery
			)
		)
	}

	# make use of the function simulateModelStochastic that returns trajectories in the correct format
	return(simulateModelStochastic(theta,state.init,times,SIR_transitions,SIR_rateFunc))

}

SIR_generateObservation <- function(model.traj, theta){

	# daily incidence needed
	incidence <- diff(model.traj$Inc)

	x <- rpois(length(incidence),lambda=theta[["rho"]]*incidence)

	model.traj$observation <- c(0,x)

	return(model.traj)
}

# function to compute log-prior
SIR_logPrior <- function(list.fitparam) {

	theta <- getParameterValues(list.fitparam)

	# check constraint, due to the constant population size, to avoid negative proportion of susceptible as initial condition
	if((theta[["pI0"]]+theta[["pR0"]]) > 1){
		return(-Inf)
	}

	return(compositeLogPrior(list.fitparam))
}

SIR_logLikelihood <- function(data, theta, model.traj){

	# daily incidence needed
	incidence <- diff(model.traj$Inc)

	# keep only data incidence corresponding to simulated times
	data <- subset(data,time%in%model.traj$time[-1]) # [-1] to remove initial simulation time

	x <- sum(dpois(x=data$Inc,lambda=theta[["rho"]]*incidence,log=TRUE))

	return(x)
}


SIR <- fitmodel(
	name="SIR",
	state.variables=c("S","I","R","Inc"),
	list.fitparam=list(R0,InfectiousPeriod,ReportingRate,proportionI0,proportionR0,PopSize),
	initialise.state=SIR_initialiseState,
	log.prior.fitparam=SIR_logPrior,
	simulate.model=SIR_simulateDeterministic,
	generate.observation=SIR_generateObservation,
	data=measles,
	log.likelihood=SIR_logLikelihood)





