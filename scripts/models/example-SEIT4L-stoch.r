seit4lStoName <-
  "stochastic SEIT4L model with daily incidence and constant population size"
seit4lStateNames <- c("S", "E", "I", "T1", "T2", "T3", "T4", "L", "Inc")

# Simulate realisation of the stochastic version of the SEIT4L model.
seit4lSimulateStochastic <- function(theta, initState, times) {
  seit4lTransitions <- list(
    c(S = -1, E = 1), # infection
    c(E = -1, I = 1, Inc = 1), # infectiousness and incidence
    c(I = -1, T1 = 1), # recovery + temporary protection
    c(T1 = -1, T2 = 1), # progression of temporary protection
    c(T2 = -1, T3 = 1), # progression of temporary protection
    c(T3 = -1, T4 = 1), # progression of temporary protection
    c(T4 = -1, L = 1), # efficient long term protection
    c(T4 = -1, S = 1) # deficient long term protection
  )

  seit4lRateFunc <- function(state, theta, t) {
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
    t3 <- state[["T3"]]
    t4 <- state[["T4"]]
    l <- state[["L"]]
    inc <- state[["Inc"]]

    n <- s + e + i + t1 + t2 + t3 + t4 + l

    return(c(
      beta * s * i / n, # new infection (S -> e)
      epsilon * e, # infectiousness and incidence (E -> i)
      nu * i, # recovery + short term protection (I -> t1)
      4 * tau * t1, # progression of temporary protection (T1 -> t2)
      4 * tau * t2, # progression of temporary protection (T2 -> t3)
      4 * tau * t3, # progression of temporary protection (T3 -> t4)
      alpha * 4 * tau * t4, # efficient long term protection (T4 -> l)
      (1 - alpha) * 4 * tau * t4 # deficient long term protection (T4 -> s)
    ))
  }

  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- fitR::simulateModelStochastic(
    theta, initState, times, seit4lTransitions, seit4lRateFunc
  )

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}

seit4lStoch <- fitmodel(
  name = seit4lStoName,
  stateNames = seit4lStateNames,
  thetaNames = seitlThetaNames,
  simulate = seit4lSimulateStochastic,
  dPrior = seitlPrior,
  rPointObs = seitlGenObsPoint,
  dPointObs = seitlPointLike
)
