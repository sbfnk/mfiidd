SEITL_deter_name <- # nolint
  "deterministic SEITL model with daily incidence and constant population size"
# note the new state Inc for the daily incidence
SEITL_stateNames <- c("S", "E", "I", "T", "L", "Inc") # nolint
SEITL_thetaNames <- c("R_0", "D_lat", "D_inf", "alpha", "D_imm", "rho") # nolint

# Solves the system of ordinary differential equations for the SEITL model.
SEITL_simulateDeterministic <- function(theta, initState, times) { # nolint
  SEITL_ode <- function(time, state, theta) { # nolint
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

    N <- S + E + I + T + L # nolint

    dS <- -beta * S * I / N + (1 - alpha) * tau * Ti
    dE <- beta * S * I / N - epsilon * E
    dI <- epsilon * E - nu * I
    dT <- nu * I - tau * Ti
    dL <- alpha * tau * Ti
    dInc <- epsilon * E

    return(list(c(dS, dE, dI, dT, dL, dInc)))
  }


  # put incidence at 0 in initState
  initState["Inc"] <- 0

  traj <- as.data.frame(deSolve::ode(
    initState, times, SEITL_ode, theta,
    method = "ode45"
  ))

  # compute incidence of each time interval
  traj$Inc <- c(0, diff(traj$Inc))

  return(traj)
}


# Generate an observed incidence under a Poisson observation process.
SEITL_genObsPoint <- function(modelPoint, theta) { # nolint
  obsPoint <- rpois(n = 1, lambda = theta[["rho"]] * modelPoint[["Inc"]])

  return(c(obs = obsPoint))
}

# Evaluate the (log of the) prior density distribution of the parameter values.
SEITL_prior <- function(theta, log = FALSE) { # nolint
  logPrior_R_0 <- dunif( # nolint
    theta[["R_0"]],
    min = 1, max = 50, log = TRUE
  )
  logPrior_latentPeriod <- dunif( # nolint
    theta[["D_lat"]],
    min = 0, max = 10, log = TRUE
  )
  logPrior_infectiousPeriod <- dunif( # nolint
    theta[["D_inf"]],
    min = 0, max = 15, log = TRUE
  )
  logPrior_temporaryImmunePeriod <- dunif( # nolint
    theta[["D_imm"]],
    min = 0, max = 50, log = TRUE
  )
  logPrior_probabilityLongTermImmunity <- dunif( # nolint
    theta[["alpha"]],
    min = 0, max = 1, log = TRUE
  )
  logPrior_reportingRate <- dunif( # nolint
    theta[["rho"]],
    min = 0, max = 1, log = TRUE
  )

  logSum <-
    logPrior_R_0 + logPrior_latentPeriod + logPrior_infectiousPeriod +
    logPrior_temporaryImmunePeriod + logPrior_probabilityLongTermImmunity +
    logPrior_reportingRate

  return(ifelse(log, logSum, exp(logSum)))
}


# Computes the (log)-likelihood of a data point given the state of the model and
# under a poisson observation process.
SEITL_pointLike <- function(dataPoint, modelPoint, theta, log = FALSE) { # nolint
  return(dpois(
    x = dataPoint[["obs"]], lambda = theta[["rho"]] * modelPoint[["Inc"]],
    log = log
  ))
}


# create fitmodel
SEITL_deter <- fitmodel( # nolint
  name = SEITL_deter_name,
  stateNames = SEITL_stateNames,
  thetaNames = SEITL_thetaNames,
  simulate = SEITL_simulateDeterministic,
  dprior = SEITL_prior,
  rPointObs = SEITL_genObsPoint,
  dPointObs = SEITL_pointLike
)
