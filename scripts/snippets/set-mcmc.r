# load the fitmodels
data(models)

# wrapper for posterior
my_posteriorTdC <- function(theta) {
  my_fitmodel <- seitlDeter
  my_initState <- c("S" = 279, "E" = 0, "I" = 2, "T" = 3, "L" = 0, "Inc" = 0)
  # note that for the SEIT4L model there are 4 state variables for the T
  # compartment
  # my_initState <- c("S" = 279, "E" = 0, "I" = 2, "T1" = 3, "T2" = 0, "T3" = 0,
  # "T4" = 0, "L" = 0, "Inc" = 0)

  return(dLogPosterior(
    fitmodel = my_fitmodel, theta = theta, initState = my_initState,
    data = fluTdc1971, margLogLike = dTrajObs, log = TRUE
  ))
}

# theta to initialise the MCMC
initTheta <- c(
  "R_0" = 2, "D_lat" = 2, "D_inf" = 2, "alpha" = 0.8, "D_imm" = 16, "rho" = 0.85
)

# diagonal elements of the covariance matrix for the Gaussian proposal
proposalSd <- c(
  "R_0" = 1, "D_lat" = 0.5, "D_inf" = 0.5, "alpha" = 0.1, "D_imm" = 2, "rho" = 0.1
)

# lower and upper limits of each parameter
lower <- c(
  "R_0" = 0, "D_lat" = 0, "D_inf" = 0, "alpha" = 0, "D_imm" = 0, "rho" = 0
)
upper <- c(
  "R_0" = Inf, "D_lat" = Inf, "D_inf" = Inf, "alpha" = 1, "D_imm" = Inf, "rho" = 1
)

# number of iterations for the MCMC
nIterations <- 5000

# additional parameters for the adaptive MCMC, see ?mcmcMh for more details
adaptSizeStart <- 100
adaptSizeCooling <- 0.999
adaptShapeStart <- 200
