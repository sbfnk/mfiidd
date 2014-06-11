# create a simple stochastic SIR model with constant population size
N <- 7e+6

# define model parameters using the fitparam class
R0 <- fitparam(name="R0",value=log(12, 10),support=c(-Inf,Inf),sd.proposal=0, prior=list(distribution="dunif",parameters=c(min=log(5, 10),max=log(25, 10))))

InfectiousPeriod <- fitparam(name="IP",value=log(2, 10))

ReportingRate <- fitparam(name="rho",value=log(1, 10),support=c(-Inf,0),sd.proposal=1e-3, prior=list(distribution="dunif",parameters=c(min=-1,max=0)))

proportionR0 <- fitparam(name="pR0",value=-0.06,support=c(-Inf,0),sd.proposal=1e-4, prior=list(distribution="dunif",parameters=c(min=log(0.5, 10),max=0)))
## proportionR0 <- fitparam(name="pR0",value=log(0.9, 10),support=c(-Inf,0),sd.proposal=5e-5, prior=list(distribution="dunif",parameters=c(min=log(0.5, 10),max=0)))
proportionR0 <- fitparam(name="pR0",value=-0.06,support=c(-Inf,0),sd.proposal=1e-4, prior=list(distribution="dunif",parameters=c(min=log(0.5, 10),max=0)))
## PopSize <- fitparam(name="N",value=log(N, 10), support=c(0, log(N, 10)), sd.proposal = 1e-3, prior=list(distribution="dunif", parameters=c(min=log(N, 10)-2, max=log(N, 10))))
PopSize <- fitparam(name="N",value=5.4, support=c(0, log(N, 10)), sd.proposal = 1e-3, prior=list(distribution="dunif", parameters=c(min=log(N, 10)-2, max=log(N, 10))))

# function to initialise the model
SIR_initialiseState <- function(theta) {

	# constant pop size
	N <- 10^theta[["N"]]

	# number of infected and immune
	I <- round(theta[["nI0"]])
	R <- round(10^theta[["pR0"]]*N)

	if(I+R>N){
		stop("Initial conditions not valid")
	}

	return(c(S=N-I-R,I=I,R=R,Inc=0))
}



SIR_simulateDeterministic <- function(theta,state.init,times) {

	SIR_ode <- function(time,state,theta) {

        ## param
		beta <- 10^theta[["R0"]]/10^theta[["IP"]]
		gamma <- 1/10^theta[["IP"]]
		N <- 10^theta[["N"]]

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
		beta <- 10^theta[["R0"]]/10^theta[["IP"]]
		nu <- 1/10^theta[["IP"]]
		N <- 10^theta[["N"]]
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
	daily.incidence <- diff(model.traj$Inc)

	x <- rpois(length(daily.incidence),lambda=10^theta[["rho"]]*daily.incidence)

	model.traj$observation <- c(0,x)

	return(model.traj)
}

# function to compute log-prior
SIR_logPrior <- function(list.fitparam) {

	theta <- getParameterValues(list.fitparam)

	# check constraint, due to the constant population size, to avoid negative proportion of susceptible as initial condition
	if((theta[["nI0"]]/10^theta[["N"]]+10^theta[["pR0"]]) > 1){
		return(-Inf)
	}

	return(compositeLogPrior(list.fitparam))
}

# TODO: use SIR mock data set
data <- data.frame(time=1,Inc=0)

SIR_logLikelihood <- function(data, theta, model.traj){

	# daily incidence needed
	daily.incidence <- diff(model.traj$Inc)

	# keep only data incidence corresponding to simulated times
	data <- subset(data,time%in%model.traj$time[-1]) # [-1] to remove initial simulation time

        x <- sum(dpois(x=data$cases,lambda=(10^theta[["rho"]]*daily.incidence),log=TRUE))

	return(x)
}

modelPosterior <- function(fitmodel) {
  return(function(theta) {
    prior <- fitmodel$log.prior(theta)

    state.init <- fitmodel$initialise.state(theta)

    trajectory <- fitmodel$simulate.model(theta, state.init, fitmodel$data$time)
    likelihood <- fitmodel$log.likelihood(fitmodel$data, theta, trajectory)

    return(list(log.dist = prior + likelihood))

  })
}

SIR <- fitmodel(
	name="SIR",
	state.variables=c("S","I","R","Inc"),
	list.fitparam=list(R0,InfectiousPeriod,ReportingRate,nI0,proportionR0,PopSize),
	initialise.state=SIR_initialiseState,
	log.prior.fitparam=SIR_logPrior,
	simulate.model=SIR_simulateDeterministic,
	generate.observation=SIR_generateObservation,
	data=measles,
	log.likelihood=SIR_logLikelihood)


SIR_posterior <- modelPosterior(SIR)

mcmc.res <- mcmcMH(target = SIR_posterior,  theta.init = SIR$theta,  gaussian.proposal = SIR$gaussian.proposal,  n.iterations = 10000)

mcmc.res$trace <- data.table(mcmc.res$trace)
mcmc.res$trace$n <- seq_len(nrow(mcmc.res$trace))

ggplot(mcmc.res$trace, aes(x = n, y = theta.pR0))+ geom_line()
ggplot(mcmc.res$trace, aes(x = n, y = theta.N))+ geom_line()
ggplot(mcmc.res$trace, aes(x = n, y = theta.rho))+ geom_line()

burn_in <- 100

theta <- SIR$theta
theta["pR0"] <- median(mcmc.res$trace[n > burn_in]$theta.pR0)
theta["N"] <- median(mcmc.res$trace[n > burn_in]$theta.N)
theta["rho"] <- median(mcmc.res$trace[n > burn_in]$theta.rho)

trajectory <- SIR$simulate.model(theta, SIR_initialiseState(theta), seq_len(nrow(measles) + 1))
trajectory$cases <- c(diff(trajectory$Inc), 0)

trajectory <- trajectory[-nrow(trajectory),  ]
trajectory$true.cases <- measles$cases * 10^theta["rho"]

saveRDS(mcmc.res, "working_mcmc.rds")
