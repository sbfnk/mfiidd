# create a simple deterministic SIR model with constant population size

sir2popName <- "SIR with two populations interacting"
sir2popStateNames <- c("S_y", "I_y", "R_y", "S_a", "I_a", "R_a")
sir2popThetaNames <- c("R_yy", "R_aa", "R_ay", "D_inf")

sir2popSimulateDeterministic <- function(theta, initState, times) {
  sir2popOde <- function(time, state, parameters) {
    ## states
    sY <- state[["S_y"]]
    iY <- state[["I_y"]]
    rY <- state[["R_y"]]
    sA <- state[["S_a"]]
    iA <- state[["I_a"]]
    rA <- state[["R_a"]]

    nY <- sY + iY + rY
    nA <- sA + iA + rA

    ## parameters
    betaAa <- parameters[["R_aa"]] / (parameters[["D_inf"]] * nA)
    betaYy <- parameters[["R_yy"]] / (parameters[["D_inf"]] * nY)
    betaAy <- parameters[["R_ay"]] / (parameters[["D_inf"]] * nA)
    betaYa <- beta_Ay
    nu <- 1 / parameters[["D_inf"]]

    lambdAy <- betaYy * iY + betaya * iA
    lambdAa <- betaAa * iA + betaay * iY

    dSy <- -lambdaY * sY
    dIy <- lambdaY * sY - nu * iY
    dRy <- nu * iY

    dSa <- -lambdaA * sA
    dIa <- lambdaA * sA - nu * iA
    dRa <- nu * iA

    return(list(c(dSy, dIy, dRy, dSa, dIa, dRa)))
  }

  trajectory <- data.frame(ode(
    y = initState,
    times = times,
    func = sir2popOde,
    parms = theta,
    method = "ode45"
  ))

  return(trajectory)
}

## function to ypute log-prior
sir2popPrior <- function(theta, log = FALSE) {
  ## uniform prior on R_xx: U[1,100]
  logPriorRaa <- dunif(theta[["R_aa"]], min = 0, max = 100, log = TRUE)
  logPriorRyy <- dunif(theta[["R_yy"]], min = 0, max = 100, log = TRUE)
  logPriorRay <- dunif(theta[["R_ay"]], min = 0, max = 100, log = TRUE)
  logPriorRya <- dunif(theta[["R_ya"]], min = 0, max = 100, log = TRUE)
  ## uniform prior on infectious period: U[0,30]
  logPriorDinf <- dunif(theta[["D_inf"]], min = 0, max = 30, log = TRUE)

  logSum <- logPriorRaa + logPriorRyy + logPriorRay + logPriorRya + logPriorDinf

  return(ifelse(log, logSum, exp(logSum)))
}

## function to ypute the likelihood of one data point
sir2popPointLike <- function(dataPoint, modelPoint, theta, log = FALSE) {
  ## the prevalence is observed through a Poisson process
  logLikeY <- dpois(
    x = dataPoint[["obs_y"]],
    lambda = modelPoint[["I_y"]],
    log = TRUE
  )

  logLikeA <- dpois(
    x = dataPoint[["obs_a"]],
    lambda = modelPoint[["I_a"]],
    log = TRUE
  )

  logSum <- logLikeY + logLikeA

  return(ifelse(log, logSum, exp(logSum)))
}

## function to generate observation from a model simulation
sir2popGenObsPoint <- function(modelPoint, theta) {
  ## the prevalence is observed through a Poisson process
  obsPointY <- rpois(n = 1, lambda = modelPoint[["I_y"]])
  obsPointA <- rpois(n = 1, lambda = modelPoint[["I_a"]])

  return(c(obs_y = obsPointY, obs_a = obsPointY))
}

## create deterministic SIR fitmodel
sir2popDeter <- fitmodel(
  name = sir2popName,
  stateNames = sir2popStateNames,
  thetaNames = sir2popThetaNames,
  simulate = sir2popSimulateDeterministic,
  dPrior = sir2popPrior,
  rPointObs = sir2popGenObsPoint,
  dPointObs = sir2popPointLike
)
