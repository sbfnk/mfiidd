# run the MCMC
my_Pmcmc <- furrr::future_map(seq_len(8), \(x) mcmcMh(
  target = my_posteriorSto,
  initTheta = initTheta,
  covmat = covmat,
  limits = list(lower = lower, upper = upper),
  nIterations = nIterations,
  adaptSizeStart = adaptSizeStart,
  adaptSizeCooling = adaptSizeCooling,
  adaptShapeStart = adaptShapeStart
), .options = furrr::furrr_options(seed = TRUE))
