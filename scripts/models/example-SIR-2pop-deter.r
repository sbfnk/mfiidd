# create a simple deterministic SIR model with constant population size

SIR_2pop_name <- "SIR with two populations interacting" # nolint
SIR_2pop_stateNames <- c("S_y", "I_y", "R_y", "S_a", "I_a", "R_a") # nolint
SIR_2pop_thetaNames <- c("R_yy", "R_aa", "R_ay", "D_inf") # nolint

SIR_2pop_simulateDeterministic <- function(theta, initState, times) { # nolint
  SIR_2pop_ode <- function(time, state, parameters) { # nolint
    ## states
    S_y <- state[["S_y"]] # nolint
    I_y <- state[["I_y"]] # nolint
    R_y <- state[["R_y"]] # nolint
    S_a <- state[["S_a"]] # nolint
    I_a <- state[["I_a"]] # nolint
    R_a <- state[["R_a"]] # nolint

    N_y <- S_y + I_y + R_y # nolint
    N_a <- S_a + I_a + R_a # nolint

    ## parameters
    beta_aa <- parameters[["R_aa"]] / (parameters[["D_inf"]] * N_a) # nolint
    beta_yy <- parameters[["R_yy"]] / (parameters[["D_inf"]] * N_y) # nolint
    beta_ay <- parameters[["R_ay"]] / (parameters[["D_inf"]] * N_a) # nolint
    beta_ya <- beta_ay # nolint
    nu <- 1 / parameters[["D_inf"]]

    lambda_y <- beta_yy * I_y + beta_ya * I_a # nolint
    lambda_a <- beta_aa * I_a + beta_ay * I_y # nolint

    dS_y <- -lambda_y * S_y # nolint
    dI_y <- lambda_y * S_y - nu * I_y # nolint
    dR_y <- nu * I_y # nolint

    dS_a <- -lambda_a * S_a # nolint
    dI_a <- lambda_a * S_a - nu * I_a # nolint
    dR_a <- nu * I_a # nolint

    return(list(c(dS_y, dI_y, dR_y, dS_a, dI_a, dR_a)))
  }

  trajectory <- data.frame(ode(
    y = initState,
    times = times,
    func = SIR_2pop_ode,
    parms = theta,
    method = "ode45"
  ))

  return(trajectory)
}

## function to ypute log-prior
SIR_2pop_prior <- function(theta, log = FALSE) { # nolint
  ## uniform prior on R_xx: U[1,100]
  logPrior_R_aa <- dunif(theta[["R_aa"]], min = 0, max = 100, log = TRUE) # nolint
  logPrior_R_yy <- dunif(theta[["R_yy"]], min = 0, max = 100, log = TRUE) # nolint
  logPrior_R_ay <- dunif(theta[["R_ay"]], min = 0, max = 100, log = TRUE) # nolint
  logPrior_R_ya <- dunif(theta[["R_ya"]], min = 0, max = 100, log = TRUE) # nolint
  ## uniform prior on infectious period: U[0,30]
  logPrior_D <- dunif(theta[["D_inf"]], min = 0, max = 30, log = TRUE) # nolint

  logSum <-
    logPrior_R_aa + logPrior_R_yy + logPrior_R_ay + logPrior_R_ya +
    logPrior_D

  return(ifelse(log, logSum, exp(logSum)))
}

## function to ypute the likelihood of one data point
SIR_2pop_pointLike <- function(data.point, modelPoint, theta, log = FALSE) { # nolint
  ## the prevalence is observed through a Poisson process
  logLike_y <- dpois( # nolint
    x = data.point[["obs_y"]],
    lambda = modelPoint[["I_y"]],
    log = TRUE
  )

  logLike_a <- dpois( # nolint
    x = data.point[["obs_a"]],
    lambda = modelPoint[["I_a"]],
    log = TRUE
  )

  logSum <- logLike_y + logLike_a

  return(ifelse(log, logSum, exp(logSum)))
}

## function to generate observation from a model simulation
SIR_2pop_genObsPoint <- function(modelPoint, theta) { # nolint
  ## the prevalence is observed through a Poisson process
  obsPoint_y <- rpois(n = 1, lambda = modelPoint[["I_y"]]) # nolint
  obsPoint_a <- rpois(n = 1, lambda = modelPoint[["I_a"]]) # nolint

  return(c(obs_y = obsPoint_y, obs_a = obsPoint_y))
}

## create deterministic SIR fitmodel
SIR_2pop_deter <- fitmodel( # nolint
  name = SIR_2pop_name,
  stateNames = SIR_2pop_stateNames,
  thetaNames = SIR_2pop_thetaNames,
  simulate = SIR_2pop_simulateDeterministic,
  dprior = SIR_2pop_prior,
  rPointObs = SIR_2pop_genObsPoint,
  dPointObs = SIR_2pop_pointLike
)
