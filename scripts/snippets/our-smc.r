## @knitr our_smc
# This is a function that takes four parameters:
# - fitmodel: a fitmodel object
# - theta: named numeric vector. Values of the parameters for which the marginal
#   log-likelihood is desired.
# - initState: named numeric vector. Initial values of the state variables.
# - data: data frame. Observation times and observed data.
# The function returns the value of the marginal log-likelihood
my_particleFilter <- function(fitmodel, theta, initState, data, nParticles) {
  ## Initialisation of the algorithm

  # Marginal log-likelihood is set to 0 and will be updated during the filtering
  # steps
  margLogLike <- 0

  # Particle states can be stored in a list
  stateParticles <- rep(list(initState), nParticles)

  # Weight: initially equal for all the particles
  # particle weight can be stored in a vector
  weightParticles <- rep(1 / nParticles, length = nParticles)

  # Initialise time variable
  currentTime <- 0

  ## Loop over observation times: resample, propagate, weight
  for (i in seq_len(nrow(data))) {
    # Extract next data point (must be a vector)
    dataPoint <- unlist(data[i, ])
    nextTime <- dataPoint["time"]

    # Resample particles according to their weights.
    # You can use the `sample` function of R
    # (normalisation of the weights is done in the function)
    indexResampled <- sample(
      x = nParticles,
      size = nParticles,
      replace = TRUE,
      prob = weightParticles
    )
    stateParticles <- stateParticles[indexResampled]

    ## Loop over particles: propagate and weight
    for (p in 1:nParticles) {
      # Extract current state of the particle
      currentStateParticle <- stateParticles[[p]]

      # Propagate the particle from current observation time
      # to the next one using the function `fitmodel$simulate`
      traj <- fitmodel$simulate(
        theta = theta,
        initState = currentStateParticle,
        times = c(currentTime, nextTime)
      )

      # Extract state of the model at next observation time
      # Also make sure that modelPoint is a vector
      modelPoint <- unlist(traj[2, fitmodel$stateNames])

      # Weight the particle with the likelihood of the observed
      # data point using the function `fitmodel$dPointObs`
      weightParticles[p] <-
        fitmodel$dPointObs(
          dataPoint = dataPoint,
          modelPoint = modelPoint,
          theta = theta
        )

      # Update state of the p particle
      stateParticles[[p]] <- modelPoint
    }

    # Increment time
    currentTime <- nextTime

    ## Increment the marginal log-likelihood
    # Add the log of the mean of the particles weights
    margLogLike <- margLogLike + log(mean(weightParticles))
  }

  ## Return marginal log-likelihood
  return(margLogLike)
}
