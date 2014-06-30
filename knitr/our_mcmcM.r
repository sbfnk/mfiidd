## @knitr our_mcmcM
## This is a function that takes four parameters:
## - target: the target distribution, a function that takes one argument
##           (a vector) and returns the (logged) value of a distribution
## - theta.init: the initial value of theta, a named vector
## - covmat.proposal: the covariance matrix of the (Gaussian) proposal distribution,
##                    in the same order as in the "target" vector
## - n.iterations: the number of iterations
## it returns an MCMC trace (value of theta and target(theta) at every MCMC step)

our_mcmcM <- function(target, theta.init, covmat.proposal, n.iterations) {

    ## initialise theta
    theta.current <- theta.init
    theta.proposed <- theta.init

    ## evaluate the function "target" at "theta.init", and assign to 
    ## a variable called target.theta.current
    target.theta.current <- target(theta.current)

    ## initialise trace data.frame
    trace <- data.frame(target = target.theta.current, t(theta.current))

    ## initialise number of accepted runs (to calculate acceptance rate later)
    accepted <- 0

    ## run MCMC for n.iteration interations
    for (i.iteration in seq_len(n.iterations)) {

        ## draw a new theta from the (multivariate Gaussian) proposal 
        ## distribution and assign to a variable called theta.proposed. 
        ## See "?rmvnorm for help.
        theta.proposed <- rmvnorm(n = 1, mean = theta.current, 
                                  sigma = covmat.proposal)

        ## evaluate the function target at the proposed theta and
        ## assign to a variable called target.theta.proposed
        target.theta.proposed <- target(theta.proposed)

        ## compute Metropolis ratio (acceptance probability). Since the
        ## multivariate Gaussian is symmetric, we don't need to consider
        ## the proposal distribution here
        log.acceptance <- target.theta.proposed - target.theta.current

        ## draw random number number between 0 and 1 using "runif" and assign to
        ## a variable called r.
        r <- runif(1)

        ## calculate acceptance ratio. This is easiest if you assume
        ## the target function to return the logarithm of the 
        ## distribution value. Assign the result to a variable called
        ## log.acceptance
        log.acceptance <- target.theta.proposed - target.theta.current

        ## test acceptance (using "exp" because we calculated the logarithm of the
        ## acceptance ratio before
        if (r < exp(log.acceptance)) {

            ## if accepted: update proposed parameter
            trace <- rbind(trace, c(target = target.theta.proposed, t(theta.proposed)))
            ## update theta
            theta.current <- theta.proposed

            ## update target
            target.theta.current <- target.theta.proposed

            ## update number of accepted proposals
            accepted <- accepted + 1
        } else {
            ## reject proposed parameter
            trace <- rbind(trace, c(target = target.theta.current, t(theta.current)))
        }

        ## print current state of chain and acceptance rate
        cat("chain:", unlist(trace[nrow(trace), ]),
            "acceptance rate:", accepted / i.iteration, "\n")

    }

    return(trace)
}
