seit2lDeterName <-
  "deterministic SEIT2L model with daily incidence and constant population size"
seit2lStateNames <- c("S", "E", "I", "T1", "T2", "L", "Inc")

seit2lSimulateDeterministic <- function(theta, initState, times) {
  seit2lOde <- function(time, state, theta) {
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

    dS <- -beta * s * i / n + (1 - alpha) * 2 * tau * t2
    dE <- beta * s * i / n - epsilon * e
    dI <- epsilon * e - nu * i
    dT1 <- nu * i - 2 * tau * t1
    dT2 <- 2 * tau * t1 - 2 * tau * t2
    dL <- alpha * 2 * tau * t2
    dInc <- epsilon * e

    return(list(c(dS, dE, dI, dT1, dT2, dL, dInc)))
  }


  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- as.data.frame(deSolve::ode(
    initState, times, seit2lOde, theta,
    method = "ode45"
  ))

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}


seit2lDeter <- fitmodel(
  name = seit2lDeterName,
  stateNames = seit2lStateNames,
  thetaNames = seitlThetaNames,
  simulate = seit2lSimulateDeterministic,
  dPrior = seitlPrior,
  rPointObs = seitlGenObsPoint,
  dPointObs = seitlPointLike
)
