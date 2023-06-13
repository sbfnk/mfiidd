# run the MCMC
my_pMcmc <- mcmcMh(
  target = my_posteriorSto,
  initTheta = initTheta,
  covmat = covmat,
  limits = list(lower = lower, upper = upper),
  nIterations = nIterations,
  adaptSizeStart = adaptSizeStart,
  adaptSizeCooling = adaptSizeCooling,
  adaptShapeStart = adaptShapeStart
)

