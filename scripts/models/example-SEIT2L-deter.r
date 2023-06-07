SEIT2L_deter_name <- # nolint
  "deterministic SEIT2L model with daily incidence and constant population size"
SEIT2L_stateNames <- c("S", "E", "I", "T1", "T2", "L", "Inc") # nolint

SEIT2L_simulateDeterministic <- function(theta, initState, times) { # nolint
  SEIT2L_ode <- function(time, state, theta) { # nolint
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

    dS <- -beta * S * I / N + (1 - alpha) * 2 * tau * T2
    dE <- beta * S * I / N - epsilon * E
    dI <- epsilon * E - nu * I
    dT1 <- nu * I - 2 * tau * T1
    dT2 <- 2 * tau * T1 - 2 * tau * T2
    dL <- alpha * 2 * tau * T2
    dInc <- epsilon * E

    return(list(c(dS, dE, dI, dT1, dT2, dL, dInc)))
  }


  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- as.data.frame(deSolve::ode(
    initState, times, SEIT2L_ode, theta,
    method = "ode45"
  ))

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}


SEIT2L_deter <- fitmodel( # nolint
  name = SEIT2L_deter_name,
  stateNames = SEIT2L_stateNames,
  thetaNames = SEITL_thetaNames,
  simulate = SEIT2L_simulateDeterministic,
  dprior = SEITL_prior,
  rPointObs = SEITL_genObsPoint,
  dPointObs = SEITL_pointLike
)
