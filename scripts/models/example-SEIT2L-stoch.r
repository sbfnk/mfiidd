SEIT2L_sto_name <- # nolint
  "stochastic SEIT2L model with daily incidence and constant population size"
SEIT2L_stateNames <- c("S", "E", "I", "T1", "T2", "L", "Inc") # nolint

# Simulate realisation of the stochastic version of the SEIT2L model.
SEIT2L_simulateStochastic <- function(theta, initState, times) { # nolint
  SEIT2L_transitions <- list( # nolint
    c(S = -1, E = 1), # infection
    c(E = -1, I = 1, Inc = 1), # infectiousness and incidence
    c(I = -1, T1 = 1), # recovery and temporary protection
    c(T1 = -1, T2 = 1), # progression of temporary protection
    c(T2 = -1, L = 1), # efficient long term protection
    c(T2 = -1, S = 1) # deficient long term protection
  )

  SEIT2L_rateFunc <- function(state, theta, t) { # nolint
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
    T1 <- state[["T1"]] # nolint
    T2 <- state[["T2"]] # nolint
    L <- state[["L"]] # nolint
    Inc <- state[["Inc"]] # nolint

    N <- S + E + I + T1 + T2 + L # nolint

    return(c(
      beta * S * I / N, # new infection (S -> E)
      epsilon * E, # infectiousness and incidence (E -> I)
      nu * I, # recovery + short term protection (I -> T1)
      2 * tau * T1, # progression of temporary protection (T1 -> T2)
      alpha * 2 * tau * T2, # efficient long term protection (T2 -> L)
      (1 - alpha) * 2 * tau * T2 # deficient long term protection (T2 -> S)
    ))
  }

  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- fitR::simulateModelStochastic(
    theta, initState, times, SEIT2L_transitions, SEIT2L_rateFunc
  )

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}


SEIT2L_stoch <- fitmodel( # nolint
  name = SEIT2L_sto_name,
  stateNames = SEIT2L_stateNames,
  thetaNames = SEITL_thetaNames,
  simulate = SEIT2L_simulateStochastic,
  dprior = SEITL_prior,
  rPointObs = SEITL_genObsPoint,
  dPointObs = SEITL_pointLike
)
