SEIT4L_deter_name <- # nolint
  "deterministic SEIT4L model with daily incidence and constant population size"
SEIT4L_stateNames <- c("S", "E", "I", "T1", "T2", "T3", "T4", "L", "Inc") # nolint

SEIT4L_simulateDeterministic <- function(theta, initState, times) { # nolint
  SEIT4L_ode <- function(time, state, theta) { # nolint
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
    T3 <- state[["T3"]] # nolint
    T4 <- state[["T4"]] # nolint
    L <- state[["L"]] # nolint
    Inc <- state[["Inc"]] # nolint

    N <- S + E + I + T1 + T2 + T3 + T4 + L # nolint

    dS <- -beta * S * I / N + (1 - alpha) * 4 * tau * T4
    dE <- beta * S * I / N - epsilon * E
    dI <- epsilon * E - nu * I
    dT1 <- nu * I - 4 * tau * T1
    dT2 <- 4 * tau * T1 - 4 * tau * T2
    dT3 <- 4 * tau * T2 - 4 * tau * T3
    dT4 <- 4 * tau * T3 - 4 * tau * T4
    dL <- alpha * 4 * tau * T4
    dInc <- epsilon * E

    return(list(c(dS, dE, dI, dT1, dT2, dT3, dT4, dL, dInc)))
  }


  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- as.data.frame(deSolve::ode(
    initState, times, SEIT4L_ode, theta,
    method = "ode45"
  ))

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}


SEIT4L_deter <- fitmodel( # nolint
  name = SEIT4L_deter_name,
  stateNames = SEIT4L_stateNames,
  thetaNames = SEITL_thetaNames,
  simulate = SEIT4L_simulateDeterministic,
  dprior = SEITL_prior,
  rPointObs = SEITL_genObsPoint,
  dPointObs = SEITL_pointLike
)
