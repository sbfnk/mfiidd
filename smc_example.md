Below, you can find an example of how to code a particle filter. Some bits are left out for you to fill in (marked "INSERT HERE"). Each "INSERT HERE" statement requires one line of code. If you struggle, you can find a link to the solution at the end of the page.


```r
# The particle filter returns an estimate of the marginal log-likelihood.
# It takes two arguments as inputs
# fitmodel: your fitmodel object
# n.particles: number of particles

my_particleFilter <- function(fitmodel, n.particles)
{

	############################################################################################
	## This function compute the marginal log.likelihood of the data (fitmodel$data) 
	## given the parameters (fitmodel$theta) using a particle filter
	############################################################################################

    # Create some temporary variables (useful to avoid repetition of long names with $)
    data <- fitmodel$data
    theta <- fitmodel$theta

    ############################################################################################
    ## Initialisation of the algorithm
    ############################################################################################

    # Marginal log-likelihood is set to 0 and will be updated during the fitering steps
    log.likelihood <- 0

    # Initialise the particles:
    initialise.particle  <-  # INSERT HERE vector of initial states returned by the function fitmodel$initialise.state
    # particle states can be stored in a list
    state.particles <- rep(list(initialise.particle),n.particles)

    weight.particles  <-  # INSERT HERE vector of weight initially equal for all the particles

    # initialise time variable
    current.time <- 0

    ############################################################################################
    # Loop over observation time: resample, propagate, weight
    ############################################################################################
    for(i in seq_len(nrow(data))){

        # extract next observtaion time from data
        next.time <- data$time[i]

        index.resampled <- # INSERT HERE resample particles indexes according to the particle weights 
        # resample states according to index
        state.particles <- state.particles[index.resampled]

        # propagate particles to the next observation time
        for(p in 1:n.particles){

            # extract current state of the particle 
            current.state.particle <- state.particles[[p]]

            # simulate from current observation time to next observation time
            traj <- # INSERT HERE

            # compute particle weight
            weight.particles[p] <- # INSERT HERE

            # Update state of the p particle
            # We take the last row of the traj data frame which corresponds to the next observation time
            # Also we make sure to only take state variables (traj contains also a time variable)
            state.particles[[p]] <- unlist(traj[2,fitmodel$state.variables])

        }

        # update time
        current.time <- next.time

        # Update marginal log-likelihood
        log.likelihood <- # INSERT HERE

    }

    # Return marginal log-likelihood
    return(log.likelihood)

}
```

You can find the solution [here](smc_solution.md)
