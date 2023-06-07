## create a simple deterministic SIR model with constant population size

SIR_name <- "SIR with constant population size" # nolint
SIR_stateNames <- c("S", "I", "R") # nolint
SIR_thetaNames <- c("R_0", "D_inf") # nolint

SIR_simulateDeterministic <- function(theta, initState, times) { # nolint
  SIR_ode <- function(time, state, parameters) { # nolint
    ## parameters
    beta <- parameters[["R_0"]] / parameters[["D_inf"]]
    nu <- 1 / parameters[["D_inf"]]

    ## states
    S <- state[["S"]] # nolint
    I <- state[["I"]] # nolint
    R <- state[["R"]] # nolint

    N <- S + I + R # nolint

    dS <- -beta * S * I / N
    dI <- beta * S * I / N - nu * I
    dR <- nu * I

    return(list(c(dS, dI, dR)))
  }

  trajectory <- data.frame(deSolve::ode(
    y = initState,
    times = times,
    func = SIR_ode,
    parms = theta,
    method = "ode45"
  ))

  return(trajectory)
}

## function to compute log-prior
SIR_prior <- function(theta, log = FALSE) { # nolint
  ## uniform prior on R_0: U[1,100]
  logPrior_R_0 <- dunif(theta[["R_0"]], min = 1, max = 100, log = TRUE) # nolint
  ## uniform prior on infectious period: U[0,30]
  logPrior_D <- dunif(theta[["D_inf"]], min = 0, max = 30, log = TRUE) # nolint

  logSum <- logPrior_R_0 + logPrior_D

  return(ifelse(log, logSum, exp(logSum)))
}

## function to compute the likelihood of one data point
SIR_pointLike <- function(dataPoint, modelPoint, theta, log = FALSE) { # nolint
  ## the prevalence is observed through a Poisson process
  return(dpois(
    x = dataPoint[["obs"]],
    lambda = modelPoint[["I"]],
    log = log
  ))
}

## function to generate observation from a model simulation
SIR_genObsPoint <- function(modelPoint, theta) { # nolint
  ## the prevalence is observed through a Poisson process
  obsPoint <- rpois(n = 1, lambda = modelPoint[["I"]])

  return(c(obs = obsPoint))
}

## create deterministic SIR fitmodel
SIR_deter <- fitmodel( # nolint
  name = SIR_name,
  stateNames = SIR_stateNames,
  thetaNames = SIR_thetaNames,
  simulate = SIR_simulateDeterministic,
  dprior = SIR_prior,
  rPointObs = SIR_genObsPoint,
  dPointObs = SIR_pointLike
)
