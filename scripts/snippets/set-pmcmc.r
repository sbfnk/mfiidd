# wrapper for posterior
my_posteriorSto <- function(theta) {
  my_fitmodel <- seit4lStoch
  my_initState <- c(
    S = 279, E = 0, I = 2, T1 = 3, T2 = 0,  T3 = 0,  T4 = 0, L = 0, Inc = 0
  )
  my_nParticles <- 128

  return(dLogPosterior(
    fitmodel = my_fitmodel,
    theta = theta,
    initState = my_initState,
    data = fluTdc1971,
    margLogLike = my_particleFilter,
    nParticles = my_nParticles
  ))
}

# load results of deterministic fit
data(mcmcTdcDeterLongRun)

# Let's use the first trace only, no need to burn or thin
trace <- mcmcSeitlInfoPriorTheta1$trace

# we will start the pMCMC at the mean posterior estimate
# of the deterministic fit
initTheta <- colMeans(trace[, seit4lStoch$thetaNames])

# and we take the empirical covariance matrix for the 
# Gaussian kernel proposal
covmat <- mcmcSeitlInfoPriorTheta1$covmatEmpirical

# lower and upper limits of each parameter
lower <- c(
  R_0 = 0, D_lat = 0 , D_inf = 0, alpha = 0, D_imm = 0, rho = 0
)
upper <- c(
  R_0 = Inf, D_lat = Inf , D_inf = Inf, alpha = 1, D_imm = Inf, rho = 1
)

# number of iterations for the MCMC
nIterations <- 50 # just a few since it takes quite a lot of time

# Here we don't adapt so that we can check the acceptance rate of the empirical covariance matrix
adaptSizeStart <- 100
adaptSizeCooling <- 0.99
adaptShapeStart <- 100

