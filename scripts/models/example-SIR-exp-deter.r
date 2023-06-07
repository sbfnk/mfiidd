# Create a simple stochastic SIR model with constant population size
# and parameters on the exponential scale
#
# This is based on the determinstic SIR model, contained in
# example-SIR-deter.r

SIR_exp_name <- paste0( # nolint
  "SIR with constant population size, parameters transformed to the ",
  "exponential scale"
)

SIR_exp_simulateDeterministic <- function(theta, initState, times) { # nolint
  SIR_ode <- function(time, state, parameters) { # nolint
    ## parameters
    beta <- exp(parameters[["R_0"]]) / exp(parameters[["D_inf"]])
    nu <- 1 / exp(parameters[["D_inf"]])

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
    y = initState, times = times, func = SIR_ode, parms = theta,
    method = "ode45"
  ))

  return(trajectory)
}

## function to compute log-prior
SIR_exp_logPrior <- function(theta, log = FALSE) { # nolint
  ## uniform prior on R_0: U[1,100]
  logPrior_R_0 <- dunif(exp(theta[["R_0"]]), min = 1, max = 100, log = TRUE) # nolint
  ## uniform prior on infectious period: U[0,30]
  logPrior_D <- dunif(exp(theta[["D_inf"]]), min = 0, max = 30, log = TRUE) # nolint

  logSum <- logPrior_R_0 + logPrior_D

  return(ifelse(log, logSum, exp(logSum)))
}

## create deterministic SIR fitmodel
SIR_exp_deter <- fitmodel( # nolint
  name = SIR_exp_name,
  stateNames = SIR_deter$stateNames,
  thetaNames = SIR_deter$thetaNames,
  simulate = SIR_exp_simulateDeterministic,
  dprior = SIR_exp_logPrior,
  rPointObs = SIR_deter$rPointObs,
  dPointObs = SIR_deter$dPointObs
)
