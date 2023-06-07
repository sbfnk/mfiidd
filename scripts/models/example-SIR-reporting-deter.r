## Create a simple stochastic SIR model with constant population size
##
## This is based on the determinstic SIR model, contained in
## example-SIR-deter.r

SIR_reporting_name <- # nolint
  "SIR with constant population size and incomplete reporting"
SIR_reporting_thetaNames <- SIR_thetaNames <- c("R_0", "D_inf", "RR") # nolint

## function to compute log-prior
SIR_logPrior <- function(theta, log = FALSE) { # nolint
  ## uniform prior on R_0: U[1,100]
  logPrior_R_0 <- dunif(theta[["R_0"]], min = 1, max = 100, log = TRUE) # nolint
  ## uniform prior on infectious period: U[0,30]
  logPrior_D_inf <- dunif(theta[["D_inf"]], min = 0, max = 30, log = TRUE) # nolint
  ## uniform prior on the reporting rate: U[0,1]
  logPrior_RR <- dunif(theta[["RR"]], min = 0, max = 1, log = TRUE) # nolint

  logSum <- logPrior_R_0 + logPrior_D_inf + logPrior_RR

  return(ifelse(log, logSum, exp(logSum)))
}

## function to compute the log-likelihood of one data point
SIR_reporting_pointLogLike <- # nolint
  function(dataPoint, modelPoint, theta, log = FALSE) {
    ## the prevalence is observed through a Poisson process with a reporting rate
    return(dpois(
      x = dataPoint[["obs"]],
      lambda = modelPoint[["I"]] * theta[["RR"]],
      log = log
    ))
  }

## function to generate observation from a model simulation
SIR_reporting_genObsPoint <- function(modelPoint, theta) { # nolint
  ## the prevalence is observed through a Poisson process
  obsPoint <- rpois(n = 1, lambda = modelPoint[["I"]] * theta[["RR"]])

  return(obsPoint)
}

## create deterministic SIR fitmodel
SIR_reporting_deter <- fitmodel( # nolint
  name = SIR_reporting_name,
  stateNames = SIR_stateNames,
  thetaNames = SIR_thetaNames,
  simulate = SIR_simulateDeterministic,
  dprior = SIR_logPrior,
  rPointObs = SIR_reporting_genObsPoint,
  dPointObs = SIR_reporting_pointLogLike
)
