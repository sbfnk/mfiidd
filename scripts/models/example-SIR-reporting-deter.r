## Create a simple stochastic SIR model with constant population size
##
## This is based on the determinstic SIR model, contained in
## example-sir-deter.r

sirReportingName <-
  "SIR with constant population size and incomplete reporting"
sirReportingThetaNames <- sirThetaNames <- c("R_0", "D_inf", "RR")

## function to compute log-prior
sirLogPrior <- function(theta, log = FALSE) {
  ## uniform prior on R_0: U[1,100]
  logPriorR0 <- dunif(theta[["R_0"]], min = 1, max = 100, log = TRUE)
  ## uniform prior on infectious period: U[0,30]
  logPriorDinf <- dunif(theta[["D_inf"]], min = 0, max = 30, log = TRUE)
  ## uniform prior on the reporting rate: U[0,1]
  logPriorRr <- dunif(theta[["RR"]], min = 0, max = 1, log = TRUE)

  logSum <- logPriorR0 + logPriorDinf + logPriorRr

  return(ifelse(log, logSum, exp(logSum)))
}

## function to compute the log-likelihood of one data point
sirReportingPointLogLike <- # nolint
  function(dataPoint, modelPoint, theta, log = FALSE) {
    ## the prevalence is observed through a Poisson process with a reporting rate
    return(dpois(
      x = dataPoint[["obs"]],
      lambda = modelPoint[["I"]] * theta[["RR"]],
      log = log
    ))
  }

## function to generate observation from a model simulation
sirReportingGenObsPoint <- function(modelPoint, theta) { # nolint
  ## the prevalence is observed through a Poisson process
  obsPoint <- rpois(n = 1, lambda = modelPoint[["I"]] * theta[["RR"]])

  return(obsPoint)
}

## create deterministic SIR fitmodel
sirReportingDeter <- fitmodel( # nolint
  name = sirReportingName,
  stateNames = sirStateNames,
  thetaNames = sirThetaNames,
  simulate = sirSimulateDeterministic,
  dPrior = sirLogPrior,
  rPointObs = sirReportingGenObsPoint,
  dPointObs = sirReportingPointLogLike
)
