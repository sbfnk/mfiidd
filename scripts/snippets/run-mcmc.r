# @knitr run_mcmc
my_mcmcTdc <- mcmcMh(
  target = my_posteriorTdc,
  initTheta = initTheta,
  proposalSd = proposalSd,
  limits = list(lower = lower,upper = upper),
  nIterations = nIterations,
  adaptSizeStart = adaptSizeStart,
  adaptSizeCooling = adaptSizeCooling,
  adaptShapeStart = adaptShapeStart
)
