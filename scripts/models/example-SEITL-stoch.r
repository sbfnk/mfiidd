SEITL_sto_name <- # nolint
  "stochastic SEITL model with daily incidence and constant population size"

# Simulate realisation of the stochastic version of the SEITL model.
SEITL_simulateStochastic <- function(theta, initState, times) { # nolint
  SEITL_transitions <- list( # nolint
    c(S = -1, E = 1), # infection
    c(E = -1, I = 1, Inc = 1), # infectiousness and incidence
    c(I = -1, T = 1), # recovery + short term protection
    c(T = -1, L = 1), # efficient long term protection
    c(T = -1, S = 1) # deficient long term protection
  )

  SEITL_rateFunc <- function(state, theta, t) { # nolint
    # param
    beta <- theta[["R_0"]] / theta[["D_inf"]]
    epsilon <- 1 / theta[["D_lat"]]
    nu <- 1 / theta[["D_inf"]]
    alpha <- theta[["alpha"]]
    tau <- 1 / theta[["D_imm"]]

    # states
    S <- state[["S"]] # nolint
    E <- state[["E"]] # nolint
    I <- state[["I"]] # nolint
    Ti <- state[["T"]] # nolint
    L <- state[["L"]] # nolint
    Inc <- state[["Inc"]] # nolint

    N <- S + E + I + Ti + L # nolint

    return(c(
      beta * S * I / N, # infection
      epsilon * E, # infectiousness and incidence
      nu * I, # recovery + short term protection
      alpha * tau * Ti, # efficient long term protection
      (1 - alpha) * tau * Ti # deficient long term protection
    ))
  }

  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- fitR::simulateModelStochastic(
    theta, initState, times, SEITL_transitions, SEITL_rateFunc
  )

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}

SEITL_stoch <- fitmodel( # nolint
  name = SEITL_sto_name,
  stateNames = SEITL_stateNames,
  thetaNames = SEITL_thetaNames,
  simulate = SEITL_simulateStochastic,
  dprior = SEITL_prior,
  rPointObs = SEITL_genObsPoint,
  dPointObs = SEITL_pointLike
)
