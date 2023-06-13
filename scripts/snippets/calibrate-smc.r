# load fitmodel, data and define initState
data(fluTdc1971)

initState <- c(
  S = 279, E = 0, I = 2, T1 = 3, T2 = 0, T3 = 0, T4 = 0, L = 0, Inc = 0
)

## pick a theta close to the mean posterior estimate of the deterministic fit
theta <- c(
  R_0 = 7, D_lat = 1, D_inf = 4, alpha = 0.5, D_imm = 10, rho = 0.65
)

## vector of number of particles to test
testNparticles <- 2**seq(0, 10)

## number of replicates
nReplicates <- 100

## vector and data frame of results
sampleLogLike <- vector("numeric", length = nReplicates)
calibrateSmc <- data.frame()

for (nParticles in testNparticles) {
  cat("Testing ", nParticles, " particles\n")
  ## start measuring time
  startTime <- Sys.time()
  for (i in seq_len(nReplicates)) {
    ## one Monte-Carlo estimate of the log-likelihood
    sampleLogLike[i] <- my_particleFilter(
      seit4lStoch, theta, initState, fluTdc1971, nParticles
    )
  }
  ## end measuring time
  endTime <- Sys.time()

  ## keep only replicate with finite log-likelihood to be able to compute the
  ## mean and sd this give us the proportion of replicates with particle
  ## depletion.
  sampleFiniteLogLike <- sampleLogLike[is.finite(sampleLogLike)]

  ans <- c(
    nParticles = nParticles,
    mean = mean(sampleFiniteLogLike),
    sd = sd(sampleFiniteLogLike),
    propDepleted = 1 - length(sampleFiniteLogLike) / length(sampleLogLike),
    days = difftime(endTime, startTime, units = "days")
  )

  calibrateSmc <- rbind(calibrateSmc, t(ans))
}
