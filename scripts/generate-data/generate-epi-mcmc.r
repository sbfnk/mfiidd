library("fitR")

data(models)
data(epi)

dLogPosteriorEpi1 <- function(theta) {
  return(dLogPosterior(
    fitmodel = sirDeter,
    theta = theta,
    initState = c(S = 999, I = 1, R = 0),
    data = epi1
  ))
}

dLogPosteriorEpi3 <- function(theta) {
  return(dLogPosterior(
    fitmodel = sirDeter,
    theta = theta,
    initState = c(S = 999, I = 1, R = 0),
    data = epi3
  ))
}

dLogPosteriorEpi4 <- function(theta) {
  return(dLogPosterior(
    fitmodel = sirReportingDeter,
    theta = theta,
    initState = c(S = 999, I = 1, R = 0),
    data = epi4
  ))
}

mcmcEpi1 <- mcmcMh(
  target = dLogPosteriorEpi1, initTheta = c(R_0 = 3, D_inf = 2),
  proposalSd = c(0.05, 0), nIterations = 10000
)
mcmcEpi3 <- mcmcMh(
  target = dLogPosteriorEpi3, initTheta = c(R_0 = 1, D_inf = 2),
  proposalSd = c(0.01, 0.1), nIterations = 1000
)
mcmcEpi4 <- mcmcMh(
  target = dLogPosteriorEpi4, initTheta = c(R_0 = 1, D_inf = 2, RR = 1),
  proposalSd = c(0.01, 0.1, 0.01), nIterations = 10000
)

save(mcmcEpi1, mcmcEpi3, mcmcEpi4, file = here::here("data", "mcmc-epi.rdata"))
