Below, you can find an example of how to code a Metropolis sampler. Some bits are left out for you to fill in (marked "INSERT HERE"). Each "INSERT HERE" statement requires one line of code. If you struggle, you can find a link to the solution below the function.


```r
## This is a function that takes four parameters:
## - target: the target distribution, a function that takes one argument
##           (a vector) and returns the (logged) value of a distribution
## - init.theta: the initial value of theta, a named vector
## - covmat.proposal: the covariance matrix of the (Gaussian) proposal distribution,
##                    in the same order as in the "target" vector
## - n.iterations: the number of iterations
## it returns an MCMC trace (value of theta and target(theta) at every MCMC step)

our_mcmcM <- function(target, init.theta, covmat.proposal, n.iterations) {

    ## initialise theta
    theta.current <- init.theta
    theta.proposed <- init.theta

    ## INSERT HERE: evaluate the function "target" at "init.theta", and assign to
    ##              a variable called target.theta.current

    ## initialise trace data.frame
    trace <- data.frame(target = target.theta.current, t(theta.current))

    ## initialise number of accepted runs (to calculate acceptance rate later)
    accepted <- 0

    ## run MCMC for n.iteration interations
    for (i.iteration in seq_len(n.iterations)) {

        ## INSERT HERE: draw a new theta from the (multivariate Gaussian) proposal
        ##              distribution and assign to a variable called theta.proposed.
        ##              See "?rmvnorm for help.

        ## INSERT HERE: evaluate the function target at the proposed theta and
        ##              assign to a variable called target.theta.proposed

        ## INSERT HERE: compute Metropolis ratio (acceptance probability).
        ##              This is easiest if you assume the target function to
        ##              return the logarithm of the distribution value. Assign
        ##              the result to a variable called log.acceptance

        ## INSERT HERE: random number number between 0 and 1 using "runif"

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
```

If you run into any problems, have a look at our [solution](mcmc_example_solution.md).
