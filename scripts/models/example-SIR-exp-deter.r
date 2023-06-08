# Create a simple stochastic SIR model with constant population size
# and parameters on the exponential scale
#
# This is based on the determinstic SIR model, contained in
# example-sir-deter.r

sirExpName <- paste0(
  "SIR with constant population size, parameters transformed to the ",
  "exponential scale"
)

sirExpSimulateDeterministic <- function(theta, initState, times) {
  sirOde <- function(time, state, parameters) {
    ## parameters
    beta <- exp(parameters[["R_0"]]) / exp(parameters[["D_inf"]])
    nu <- 1 / exp(parameters[["D_inf"]])

    ## states
    s <- state[["S"]]
    i <- state[["I"]]
    r <- state[["R"]]

    n <- s + i + r # nolint

    dS <- -beta * s * i / n
    dI <- beta * s * i / n - nu * i
    dR <- nu * i

    return(list(c(dS, dI, dR)))
  }

  trajectory <- data.frame(deSolve::ode(
    y = initState, times = times, func = sirOde, parms = theta,
    method = "ode45"
  ))

  return(trajectory)
}

## function to compute log-prior
sirExpLogPrior <- function(theta, log = FALSE) { # nolint
  ## uniform prior on R_0: U[1,100]
  logPriorR0 <- dunif(exp(theta[["R_0"]]), min = 1, max = 100, log = TRUE)
  ## uniform prior on infectious period: U[0,30]
  logPriorDinf <- dunif(exp(theta[["D_inf"]]), min = 0, max = 30, log = TRUE)

  logSum <- logPriorR0 + logPriorDinf

  return(ifelse(log, logSum, exp(logSum)))
}

## create deterministic SIR fitmodel
sirExpDeter <- fitmodel( # nolint
  name = sirExpName,
  stateNames = sirDeter$stateNames,
  thetaNames = sirDeter$thetaNames,
  simulate = sirExpSimulateDeterministic,
  dPrior = sirExpLogPrior,
  rPointObs = sirDeter$rPointObs,
  dPointObs = sirDeter$dPointObs
)
