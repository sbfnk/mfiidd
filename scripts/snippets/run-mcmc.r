# run the MCMC
my_mcmcTdC <- mcmcMh(
  target = my_posteriorTdC,
  initTheta = initTheta,
  proposalSd = proposalSd,
  limits = list(lower = lower,upper = upper),
  nIterations = nIterations,
  adaptSizeStart = adaptSizeStart,
  adaptSizeCooling = adaptSizeCooling,
  adaptShapeStart = adaptShapeStart
)
