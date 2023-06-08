# Create a simple stochastic SIR model with constant population size
#
# This is based on the determinstic SIR model, contained in
# example-sir-deter.r

sirStochasticName <- "stochastic SIR with constant population size"

sirSimulateStochastic <- function(theta, initState, times) {
  ## transitions
  sirTransitions <- list(
    c(S = -1, I = 1), # infection
    c(I = -1, R = 1) # recovery
  )

  ## rates
  sirRateFunc <- function(x, parameters, t) {
    beta <- parameters[["R_0"]] / parameters[["D_inf"]]
    nu <- 1 / parameters[["D_inf"]]

    s <- x[["S"]]
    i <- x[["I"]]
    r <- x[["R"]]

    n <- s + i + r

    return(c(
      beta * s * i / n, # infection
      nu * i # recovery
    ))
  }

  # make use of the function simulateModelStochastic that
  # returns trajectories in the correct format
  return(fitR::simulateModelStochastic(
    theta, initState, times, sirTransitions, sirRateFunc
  ))
}

# create stochastic SIR fitmodel
sirStoch <- fitmodel( # nolint
  name = sirStochasticName,
  stateNames = sirStateNames,
  thetaNames = sirThetaNames,
  simulate = sirSimulateStochastic,
  dPrior = sirPrior,
  rPointObs = sirGenObsPoint,
  dPointObs = sirPointLike
)
