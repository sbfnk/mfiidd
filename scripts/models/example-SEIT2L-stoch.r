seit2lStoName <-
  "stochastic SEIT2L model with daily incidence and constant population size"
seit2lStateNames <- c("S", "E", "I", "T1", "T2", "L", "Inc")

# Simulate realisation of the stochastic version of the SEIT2L model.
seit2lSimulateStochastic <- function(theta, initState, times) {
  seit2lTransitions <- list(
    c(S = -1, E = 1), # infection
    c(E = -1, I = 1, Inc = 1), # infectiousness and incidence
    c(I = -1, T1 = 1), # recovery and temporary protection
    c(T1 = -1, T2 = 1), # progression of temporary protection
    c(T2 = -1, L = 1), # efficient long term protection
    c(T2 = -1, S = 1) # deficient long term protection
  )

  seit2lRateFunc <- function(state, theta, t) {
    # param
    beta <- theta[["R_0"]] / theta[["D_inf"]]
    epsilon <- 1 / theta[["D_lat"]]
    nu <- 1 / theta[["D_inf"]]
    alpha <- theta[["alpha"]]
    tau <- 1 / theta[["D_imm"]]

    # states
    s <- state[["S"]]
    e <- state[["E"]]
    i <- state[["I"]]
    t1 <- state[["T1"]]
    t2 <- state[["T2"]]
    l <- state[["L"]]
    inc <- state[["Inc"]]

    n <- s + e + i + t1 + t2 + l

    return(c(
      beta * s * i / n, # new infection (S -> e)
      epsilon * e, # infectiousness and incidence (E -> i)
      nu * i, # recovery + short term protection (I -> t1)
      2 * tau * t1, # progression of temporary protection (T1 -> t2)
      alpha * 2 * tau * t2, # efficient long term protection (T2 -> l)
      (1 - alpha) * 2 * tau * t2 # deficient long term protection (T2 -> s)
    ))
  }

  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- fitR::simulateModelStochastic(
    theta, initState, times, seit2lTransitions, seit2lRateFunc
  )

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}


seit2lStoch <- fitmodel(
  name = seit2lStoName,
  stateNames = seit2lStateNames,
  thetaNames = seitlThetaNames,
  simulate = seit2lSimulateStochastic,
  dPrior = seitlPrior,
  rPointObs = seitlGenObsPoint,
  dPointObs = seitlPointLike
)
