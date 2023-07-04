# @knitr epi3_mcmc
mcmcEpi3 <- my_mcmcMh(
  target = my_dLogPosteriorEpi3,
  initTheta = c(R_0 = 1, D_inf = 2),
  proposalSd = c(0.01, 0.1),
  nIterations = 10000
)
