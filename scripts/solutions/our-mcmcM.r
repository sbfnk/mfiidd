## @knitr our_mcmcM
## This is a function that takes four parameters:
## - target: the target distribution, a function that takes one argument
##           (a vector) and returns the (logged) value of a distribution
## - initTheta: the initial value of theta, a named vector
## - covmatProposal: the covariance matrix of the (Gaussian) proposal
##                    distribution, in the same order as in the "target" vector
## - nIterations: the number of iterations
## it returns an MCMC trace (value of theta and target(theta) at every MCMC
## step)

our_mcmcM <- function(target, initTheta, covmatProposal, nIterations) { # nolint
  ## initialise theta
  thetaCurrent <- initTheta
  thetaProposed <- initTheta

  ## evaluate the function "target" at "initTheta", and assign to
  ## a variable called targetThetaCurrent
  targetThetaCurrent <- target(thetaCurrent)

  ## initialise trace data.frame
  trace <- data.frame(target = targetThetaCurrent, t(thetaCurrent))

  ## initialise number of accepted runs (to calculate acceptance rate later)
  accepted <- 0

  ## run MCMC for nIteration interations
  for (iIteration in seq_len(nIterations)) {
    ## draw a new theta from the (multivariate Gaussian) proposal
    ## distribution and assign to a variable called thetaProposed.
    ## See "?rmvnorm for help.
    thetaProposed <- tmvtnorm::rmvnorm(
      n = 1, mean = thetaCurrent,
      sigma = covmatProposal
    )

    ## evaluate the function target at the proposed theta and
    ## assign to a variable called targetThetaProposed
    targetThetaProposed <- target(thetaProposed)

    ## compute Metropolis ratio (acceptance probability). Since the
    ## multivariate Gaussian is symmetric, we don't need to consider
    ## the proposal distribution here
    logAcceptance <- targetThetaProposed - targetThetaCurrent

    ## draw random number number between 0 and 1 using "runif" and assign to
    ## a variable called r.
    r <- runif(1)

    ## calculate acceptance ratio. This is easiest if you assume
    ## the target function to return the logarithm of the
    ## distribution value. Assign the result to a variable called
    ## logAcceptance
    logAcceptance <- targetThetaProposed - targetThetaCurrent

    ## test acceptance (using "exp" because we calculated the logarithm of the
    ## acceptance ratio before
    if (r < exp(logAcceptance)) {
      ## if accepted: update proposed parameter
      trace <- rbind(trace, c(target = targetThetaProposed, t(thetaProposed)))
      ## update theta
      thetaCurrent <- thetaProposed

      ## update target
      targetThetaCurrent <- targetThetaProposed

      ## update number of accepted proposals
      accepted <- accepted + 1
    } else {
      ## reject proposed parameter
      trace <- rbind(trace, c(target = targetThetaCurrent, t(thetaCurrent)))
    }

    ## print current state of chain and acceptance rate
    cat(
      "chain:", unlist(trace[nrow(trace), ]),
      "acceptance rate:", accepted / iIteration, "\n"
    )
  }

  return(trace)
}
