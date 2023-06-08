# @knitr generate_samples
trace <- my_mcmcMh(
  target = my_dLogPosteriorR0Epi1, # target distribution
  initTheta = 1, # intial parameter guess
  proposalSd = 0.1, # standard deviation of
  # Gaussian proposal: 0.1
  nIterations = 1000
) # number of iterations
