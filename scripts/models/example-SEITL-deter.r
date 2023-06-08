seitlDeterName <-
  "deterministic SEITL model with daily incidence and constant population size"
# note the new state Inc for the daily incidence
seitlStateNames <- c("S", "E", "I", "T", "L", "Inc")
seitlThetaNames <- c("R_0", "D_lat", "D_inf", "alpha", "D_imm", "rho")

# Solves the system of ordinary differential equations for the SEITL model.
seitlSimulateDeterministic <- function(theta, initState, times) {
  seitlOde <- function(time, state, theta) {
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

    n <- s + e + i + t + l

    dS <- -beta * s * i / n + (1 - alpha) * tau * t
    dE <- beta * s * i / n - epsilon * e
    dI <- epsilon * e - nu * i
    dT <- nu * i - tau * ti
    dL <- alpha * tau * ti
    dInc <- epsilon * e

    return(list(c(dS, dE, dI, dT, dL, dInc)))
  }


  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- as.data.frame(deSolve::ode(
    initState, times, seitlOde, theta,
    method = "ode45"
  ))

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}


# Generate an observed incidence under a Poisson observation process.
seitlGenObsPoint <- function(modelPoint, theta) {
  obsPoint <- rpois(n = 1, lambda = theta[["rho"]] * modelPoint[["Inc"]])

  return(c(obs = obsPoint))
}

# Evaluate the (log of the) prior density distribution of the parameter values.
seitlPrior <- function(theta, log = FALSE) {
  logPriorR0 <- dunif(
    theta[["R_0"]],
    min = 1, max = 50, log = TRUE
  )
  logPriorLatentPeriod <- dunif(
    theta[["D_lat"]],
    min = 0, max = 10, log = TRUE
  )
  logPriorInfectiousPeriod <- dunif(
    theta[["D_inf"]],
    min = 0, max = 15, log = TRUE
  )
  logPriorTemporaryImmunePeriod <- dunif(
    theta[["D_imm"]],
    min = 0, max = 50, log = TRUE
  )
  logPriorProbabilityLongTermImmunity <- dunif(
    theta[["alpha"]],
    min = 0, max = 1, log = TRUE
  )
  logPriorReportingRate <- dunif(
    theta[["rho"]],
    min = 0, max = 1, log = TRUE
  )

  logSum <-
    logPriorR0 + logPriorLatentPeriod + logPriorInfectiousPeriod +
    logPriorTemporaryImmunePeriod + logPriorProbabilityLongTermImmunity +
    logPriorReportingRate

  return(ifelse(log, logSum, exp(logSum)))
}


# Computes the (log)-likelihood of a data point given the state of the model and
# under a poisson observation process.
seitlPointLike <- function(dataPoint, modelPoint, theta, log = FALSE) {
  return(dpois(
    x = dataPoint[["obs"]], lambda = theta[["rho"]] * modelPoint[["Inc"]],
    log = log
  ))
}


# create fitmodel
seitlDeter <- fitmodel(
  name = seitlDeterName,
  stateNames = seitlStateNames,
  thetaNames = seitlThetaNames,
  simulate = seitlSimulateDeterministic,
  dPrior = seitlPrior,
  rPointObs = seitlGenObsPoint,
  dPointObs = seitlPointLike
)
