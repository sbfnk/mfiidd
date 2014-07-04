

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
    # initial state: vector calculated by the function fitmodel$initialise.state evaluated at theta
    initialise.particle  <- fitmodel$initialise.state(theta)

    # particle states can be stored in a list
    state.particles <- rep(list(initialise.particle),n.particles)
    # weight: initially equal for all the particles 
    # particle weight can be stored in a vector
    weight.particles <- rep(1/n.particles,length=n.particles)

    # initialise time variable
    current.time <- 0

    ############################################################################################
    # Loop over observation time: resample, propagate, weight
    ############################################################################################
    for(i in seq_len(nrow(data))){

        # extract next observtaion time from data
        next.time <- data$time[i]

        # resample particles according to their weight (normalization of the weight is done in the function sample())
        index.resampled <- sample(x=n.particles,size=n.particles,replace=T,prob=weight.particles)
        state.particles <- state.particles[index.resampled]

        # propagate particles to the next observation time
        for(p in 1:n.particles){
            
            # extract current state of the particle 
            current.state.particle <- state.particles[[p]]

            # simulate from current observation time to next observation time
            traj <- fitmodel$simulate.model(theta=theta,state.init=current.state.particle,times=c(current.time,next.time))

            # compute particle weight
            weight.particles[p] <- exp(fitmodel$log.likelihood( data=data, model.traj=traj, theta= theta))

            # Update state of the p particle
            # You need to take last row of the traj data frame which corresponds to the next observation time
            # Also make sure to only take state variables (traj contains also a time variable)
            state.particles[[p]] <- unlist(traj[2,fitmodel$state.variables])

        }

        # update time
        current.time <- next.time

        # Update marginal log-likelihood
        log.likelihood <- log.likelihood + log(mean(weight.particles))

    }

    # Return marginal log-likelihood
    return(log.likelihood)

}
```