seitlStochName <-
  "stochastic SEITL model with daily incidence and constant population size"

# Simulate realisation of the stochastic version of the SEITL model.
seitlSimulateStochastic <- function(theta, initState, times) {
  seitlTransitions <- list(
    c(S = -1, E = 1), # infection
    c(E = -1, I = 1, Inc = 1), # infectiousness and incidence
    c(I = -1, T = 1), # recovery + short term protection
    c(T = -1, L = 1), # efficient long term protection
    c(T = -1, S = 1) # deficient long term protection
  )

  seitlRateFunc <- function(state, theta, t) {
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
    t <- state[["T"]]
    l <- state[["L"]]
    inc <- state[["Inc"]]

    n <- s + e + i +t + l

    return(c(
      beta * s * i / n, # infection
      epsilon * e, # infectiousness and incidence
      nu * i, # recovery + short term protection
      alpha * tau * t, # efficient long term protection
      (1 - alpha) * tau * t # deficient long term protection
    ))
  }

  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- fitR::simulateModelStochastic(
    theta, initState, times, seitlTransitions, seitlRateFunc
  )

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}

seitlStoch <- fitmodel(
  name = seitlStochName,
  stateNames = seitlStateNames,
  thetaNames = seitlThetaNames,
  simulate = seitlSimulateStochastic,
  dPrior = seitlPrior,
  rPointObs = seitlGenObsPoint,
  dPointObs = seitlPointLike
)
