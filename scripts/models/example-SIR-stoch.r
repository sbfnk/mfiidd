# Create a simple stochastic SIR model with constant population size
#
# This is based on the determinstic SIR model, contained in
# example-SIR-deter.r

SIR_stochastic_name <- "stochastic SIR with constant population size" # nolint

SIR_simulateStochastic <- function(theta, initState, times) { # nolint
  ## transitions
  SIR_transitions <- list( # nolint
    c(S = -1, I = 1), # infection
    c(I = -1, R = 1) # recovery
  )

  ## rates
  SIR_rateFunc <- function(x, parameters, t) { # nolint
    beta <- parameters[["R_0"]] / parameters[["D_inf"]]
    nu <- 1 / parameters[["D_inf"]]

    S <- x[["S"]] # nolint
    I <- x[["I"]] # nolint
    R <- x[["R"]] # nolint

    N <- S + I + R # nolint

    return(c(
      beta * S * I / N, # infection
      nu * I # recovery
    ))
  }

  # make use of the function simulateModelStochastic that
  # returns trajectories in the correct format
  return(fitR::simulateModelStochastic(
    theta, initState, times, SIR_transitions, SIR_rateFunc
  ))
}

# create stochastic SIR fitmodel
SIR_stoch <- fitmodel( # nolint
  name = SIR_stochastic_name,
  stateNames = SIR_stateNames,
  thetaNames = SIR_thetaNames,
  simulate = SIR_simulateStochastic,
  dprior = SIR_prior,
  rPointObs = SIR_genObsPoint,
  dPointObs = SIR_pointLike
)
