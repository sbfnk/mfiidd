# @knitr our_mcmcMh
# This is a function that takes four parameters:
# - target: the target distribution, a function that takes one
#   argument (a number) and returns the (logged) value of a
#   distribution
# - initTheta: the initial value of theta, a number
# - proposalSd: the standard deviation of (Gaussian) proposal
#   distribution
# - nIterations: the number of iterations
# The function returns a vector of samples of theta from the target
# distribution
my_mcmcMh <- function(target, initTheta, proposalSd, nIterations) {
  # evaluate the function "target" at "initTheta", and assign to
  # a variable called targetThetaCurrent.
  targetThetaCurrent <- target(initTheta)

  # initialise variables to store the current value of theta, the
  # vector of samples, and the number of accepted runs
  thetaCurrent <- initTheta
  samples <- thetaCurrent
  accepted <- 0

  # run MCMC for nIteration interations
  for (iIteration in seq_len(nIterations)) {
    # draw a new theta from the (Gaussian) proposal distribution
    # and assign to a variable called thetaProposed.
    # See "?rnorm for more information
    # Note that this step is vectorized for any arbitratry theta
    # which will be useful when we will sample from a multivariate
    # target distribution
    thetaProposed <- rnorm(
      n = length(thetaCurrent),
      mean = thetaCurrent,
      sd = proposalSd
    )

    # Note that 'rnorm' returns an unnamed vector, but the functions of
    # 'fitmodel' need a named parameter vector. We therefore set
    # the names of thetaProposed to be the same as the names of
    # thetaCurrent
    names(thetaProposed) <- names(thetaCurrent)

    # evaluate the function target at the proposed theta and
    # assign to a variable called targetThetaProposed
    targetThetaProposed <- target(thetaProposed)

    # compute Metropolis-Hastings ratio (acceptance probability). Since
    # the multivariate Gaussian is symmetric, we don't need to consider
    # the proposal distribution here
    logAcceptance <- targetThetaProposed - targetThetaCurrent

    # draw random number number between 0 and 1 using "runif" and assign to
    # a variable called r.
    r <- runif(1)

    # test acceptance by comparing the random number to the
    # Metropolis-Hastings ratio (acceptance probability) (using
    # "exp" because we calculated the logarithm of the
    # Metropolis-Hastings ratio before)
    if (r < exp(logAcceptance)) {
      # if accepted:
      # change the current value of theta to the proposed theta
      thetaCurrent <- thetaProposed

      # updated the current value of the target
      targetThetaCurrent <- targetThetaProposed

      # update number of accepted proposals
      accepted <- accepted + 1
    }

    # add the current theta to the vector of samples
    # Note that we use `rbind` in order to deal with multivariate
    # target. So if `theta` is a vector then `samples` is a matrix.
    samples <- rbind(samples, thetaCurrent, deparse.level = 0)

    # print current state of chain and acceptance rate
    # use paste() to deal with the case where `theta` is a vector
    message(
      "iteration: ", iIteration, ", chain:", paste(thetaCurrent, collapse = " "),
      ", acceptance rate:", accepted / iIteration
    )
  }

  # return the trace of the chain (i.e., the vector of samples)
  return(samples)
}
