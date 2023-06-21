library("fitR")

data(models)
data(epi)

dLogPosteriorEpi1 <- function(theta) {
  return(dLogPosterior(
    fitmodel = sirDeter,
    theta = theta,
    initState = c(S = 999, I = 1, R = 0),
    data = epi1
  )$logDensity)
}

dLogPosteriorEpi3 <- function(theta) {
  return(dLogPosterior(
    fitmodel = sirDeter,
    theta = theta,
    initState = c(S = 999, I = 1, R = 0),
    data = epi3
  )$logDensity)
}

dLogPosteriorEpi4 <- function(theta) {
  return(dLogPosterior(
    fitmodel = sirReportingDeter,
    theta = theta,
    initState = c(S = 999, I = 1, R = 0),
    data = epi4
  )$logDensity)
}

source(here::here("scripts", "snippets", "our-mcmcMh.r"))

mcmcEpi1 <- my_mcmcMh(
  target = dLogPosteriorEpi1, initTheta = c(R_0 = 3, D_inf = 2),
  proposalSd = c(0.05, 0), nIterations = 10000
)

my_dLogPosteriorEpi3 <- dLogPosteriorEpi3
source(here::here("scripts", "snippets", "epi3-mcmc.r"))

mcmcEpi4 <- my_mcmcMh(
  target = dLogPosteriorEpi4, initTheta = c(R_0 = 1, D_inf = 2, RR = 1),
  proposalSd = c(0.01, 0.1, 0.01), nIterations = 10000
)

save(mcmcEpi1, mcmcEpi3, mcmcEpi4, file = here::here("data", "mcmcEpi.rdata"))
