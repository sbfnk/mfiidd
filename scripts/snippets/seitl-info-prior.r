seitlInfoPrior <- function(theta, log = FALSE) {
  logPriorR0 <- dunif(
    theta[["R_0"]],
    min = 1, max = 50, log = TRUE
  )
  logPriorLatentPeriod <- dnorm(
    theta[["D_lat"]],
    mean = 2, sd = 1, log = TRUE
  )
  logPriorInfectiousPeriod <- dnorm(
    theta[["D_inf"]],
    mean = 2, sd = 1, log = TRUE
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

