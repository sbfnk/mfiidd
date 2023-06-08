## create a simple deterministic SIR model with constant population size

sirName <- "SIR with constant population size"
sirStateNames <- c("S", "I", "R")
sirThetaNames <- c("R_0", "D_inf")

sirSimulateDeterministic <- function(theta, initState, times) {
  sirOde <- function(time, state, parameters) {
    ## parameters
    beta <- parameters[["R_0"]] / parameters[["D_inf"]]
    nu <- 1 / parameters[["D_inf"]]

    ## states
    s <- state[["S"]]
    i <- state[["I"]]
    r <- state[["R"]]

    n <- s + i + r

    dS <- -beta * s * i / n
    dI <- beta * s * i / n - nu * i
    dR <- nu * i

    return(list(c(dS, dI, dR)))
  }

  trajectory <- data.frame(deSolve::ode(
    y = initState,
    times = times,
    func = sirOde,
    parms = theta,
    method = "ode45"
  ))

  return(trajectory)
}

## function to compute log-prior
sirPrior <- function(theta, log = FALSE) {
  ## uniform prior on R_0: U[1,100]
  logPriorR0 <- dunif(theta[["R_0"]], min = 1, max = 100, log = TRUE)
  ## uniform prior on infectious period: U[0,30]
  logPriorDinf <- dunif(theta[["D_inf"]], min = 0, max = 30, log = TRUE)

  logSum <- logPriorR0 + logPriorDinf

  return(ifelse(log, logSum, exp(logSum)))
}

## function to compute the likelihood of one data point
sirPointLike <- function(dataPoint, modelPoint, theta, log = FALSE) {
  ## the prevalence is observed through a Poisson process
  return(dpois(
    x = dataPoint[["obs"]],
    lambda = modelPoint[["I"]],
    log = log
  ))
}

## function to generate observation from a model simulation
sirGenObsPoint <- function(modelPoint, theta) {
  ## the prevalence is observed through a Poisson process
  obsPoint <- rpois(n = 1, lambda = modelPoint[["I"]])

  return(c(obs = obsPoint))
}

## create deterministic SIR fitmodel
sirDeter <- fitmodel(
  name = sirName,
  stateNames = sirStateNames,
  thetaNames = sirThetaNames,
  simulate = sirSimulateDeterministic,
  dPrior = sirPrior,
  rPointObs = sirGenObsPoint,
  dPointObs = sirPointLike
)
