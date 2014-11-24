## @knitr our_smc
# This is a function that takes four parameters:
# - fitmodel: a fitmodel object
# - theta: named numeric vector. Values of the parameters for which the marginal log-likelihood is desired.
# - init.state: named numeric vector. Initial values of the state variables.
# - data: data frame. Observation times and observed data.
# The function returns the value of the marginal log-likelihood
my_particleFilter <- function(fitmodel, theta, init.state, data, n.particles) {

    ## Initialisation of the algorithm

    # Marginal log-likelihood is set to 0 and will be updated during the filtering steps
    margLogLike <- 0

    # Particle states can be stored in a list
    state.particles  <- rep(list(init.state),n.particles)

    # Weight: initially equal for all the particles 
    # particle weight can be stored in a vector
    weight.particles <- rep(1/n.particles,length=n.particles)

    # Initialise time variable
    current.time <- 0

    ## Loop over observation times: resample, propagate, weight
    for(i in seq_len(nrow(data))){

        # Extract next data point (must be a vector)
        data.point <- unlist(data[i, ])
        next.time <- data.point["time"]

        # Resample particles according to their weights. 
        # You can use the `sample` function of R (normalization of the weights is done in the function)
        index.resampled <- sample(x=n.particles,size=n.particles,replace=TRUE,prob=weight.particles)
        state.particles <- state.particles[index.resampled]

        ## Loop over particles: propagate and weight
        for(p in 1:n.particles){

            # Extract current state of the particle 
            current.state.particle <- state.particles[[p]]

            # Propagate the particle from current observation time 
            # to the next one using the function `fitmodel$simulate`
            traj <- fitmodel$simulate(theta=theta,init.state=current.state.particle,times=c(current.time,next.time))

            # Extract state of the model at next observation time
            # Also make sure that model.point is a vector
            model.point <- unlist(traj[2,fitmodel$state.names])

            # Weight the particle with the likelihood of the observed 
            # data point using the function `fitmodel$pointLogLike`
            weight.particles[p] <- exp(fitmodel$pointLogLike(data.point=data.point, model.point=model.point, theta=theta))

            # Update state of the p particle
            state.particles[[p]] <- model.point

        }

        # Increment time
        current.time <- next.time

        ## Increment the marginal log-likelihood
        # Add the log of the mean of the particles weights
        margLogLike <- margLogLike + log(mean(weight.particles))
    }

    ## Return marginal log-likelihood
    return(margLogLike)

}