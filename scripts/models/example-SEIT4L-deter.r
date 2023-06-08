seit4lDeterName <-
  "deterministic SEIT4L model with daily incidence and constant population size"
seit4lStateNames <- c("S", "E", "I", "T1", "T2", "T3", "T4", "L", "Inc")

seit4lSimulateDeterministic <- function(theta, initState, times) {
  seit4lOde <- function(time, state, theta) {
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

    dS <- -beta * s * i / n + (1 - alpha) * 4 * tau * t4
    dE <- beta * s * i / n - epsilon * e
    dI <- epsilon * e - nu * i
    dT1 <- nu * i - 4 * tau * t1
    dT2 <- 4 * tau * t1 - 4 * tau * t2
    dT3 <- 4 * tau * t2 - 4 * tau * t3
    dT4 <- 4 * tau * t3 - 4 * tau * t4
    dL <- alpha * 4 * tau * t4
    dInc <- epsilon * e

    return(list(c(dS, dE, dI, dT1, dT2, dT3, dT4, dL, dInc)))
  }


  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- as.data.frame(deSolve::ode(
    initState, times, seit4lOde, theta,
    method = "ode45"
  ))

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}


seit4lDeter <- fitmodel(
  name = seit4lDeterName,
  stateNames = seit4lStateNames,
  thetaNames = seitlThetaNames,
  simulate = seit4lSimulateDeterministic,
  dPrior = seitlPrior,
  rPointObs = seitlGenObsPoint,
  dPointObs = seitlPointLike
)
