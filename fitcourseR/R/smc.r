#'Bootstrap Particle Filter for fitmodel object
#'
#'The bootstrap particle filter returns an estimate of the marginal log-likelihood \eqn{L = p(y(t_{1:T})|\theta)}
#'as well as the set of filtered trajectories and their respective weights at the last observation time \eqn{\omega(t_T)=p(y(t_T)|\theta)}.
#' @param n.particles number of particles
#' @param progress if \code{TRUE} progression of the filter is displayed in the console.
#' @param n.cores number of cores on which propogation of the particles is parallelised. By default no parallelisation (\code{n.cores=1}). If \code{NULL}, set to the value returned by \code{\link[parallel]{detectCores}}.
#' @inheritParams marginalLogLikelihoodDeterministic
#' @note An unbiased state sample \eqn{x(t_{0:T}) ~ p(X(t_{0:T})|\theta,y(t_{0:T}))} can be obtained by sampling the set of trajectories \code{traj} with probability \code{traj.weight}.
#' @export
#' @seealso plotSMC
#' @import parallel doParallel
#' @return A list of 3 elements:
#' \itemize{
#' \item \code{log.likelihood} the marginal log-likelihood of the theta.
#' \item \code{traj} a list of size \code{n.particles} with all filtered trajectories.
#' \item \code{traj.weight} a vector of size \code{n.particles} with the normalised weight of the filtered trajectories.
#' }
bootstrapParticleFilter <- function(fitmodel, n.particles, progress = FALSE, n.cores = 1)
{

    if(is.null(n.cores)){
        n.cores <- detectCores()
        # cat("SMC runs on ",n.cores," cores\n")
    }

    if(n.cores > 1){
        registerDoParallel(cores=n.cores)
    }

    ## compute the log.likelihood using a particle filter

    # useful variable (avoid repetition of long names)
    data <- fitmodel$data
    theta <- fitmodel$theta

    # initialisation

    # marginal log-likelihood of the theta
    log.likelihood <- 0

    # initial state of particles
    initialise.particle  <- fitmodel$initialise.state(theta)
    current.state.particles <- rep(list(initialise.particle),n.particles)

    # filtered trajectories (just add time variable to initial state)
    traj.particles <- rep(list(data.frame(t(c(time=0,initialise.particle)))),n.particles)

    # weight of particles
    weight.particles <- rep(1/n.particles,length=n.particles)

    if(progress){
        # help to visualise progression of the filter
        progress.bar <- txtProgressBar(min=1, max= nrow(data))
    }

    # particle filter
    for(i in seq_len(nrow(data))){

        # initial + observation times
        times <- c(ifelse(i==1,0,data$time[i-1]),data$time[i])

        if(!all(weight.particles==0)){
            # resample particles according to their weight (normalization is done in the function sample())
            index.resampled <- sample(x=n.particles,size=n.particles,replace=T,prob=weight.particles)
        }else{
            warning("All particles depleted at step ",i," of SMC. Return log.likelihood = -Inf for theta set: ",paste(getParameterValues(theta),collapse=", "))
            return(list(log.likelihood=-Inf,traj=NA,traj.weight=NA))
        }

        # update traj and current state after resampling
        traj.particles <- traj.particles[index.resampled]
        current.state.particles <- current.state.particles[index.resampled]

        # propagate particles (this for loop could be parallelized)
        propagate <- llply(current.state.particles,function(current.state) {

            # simulate from previous observation to current observation time
            traj <- fitmodel$simulate.model(theta=theta,state.init=unlist(current.state),times=times)

            # compute particle weight
            weight <- exp(fitmodel$log.likelihood( data=data, model.traj=traj, theta= theta))

            return(list(traj=traj[-1,],weight=weight))

        },.parallel=(n.cores > 1))

        # collect parallel jobs
        current.state.particles <- llply(propagate,function(x) {x$traj[fitmodel$state.variables]})
        weight.particles <- unlist(llply(propagate,function(x) {x$weight}))
        traj.particles <- llply(seq_along(propagate),function(j) {rbind(traj.particles[[j]],propagate[[j]]$traj)})

        # update marginal log-likelihood
        log.likelihood <- log.likelihood + log(mean(weight.particles))

        if(progress){
            # advance progress bar
            setTxtProgressBar(progress.bar, i)
        }
    }

    if(progress){
        close(progress.bar)
    }

    # return marginal log-likelihood, filtered trajectories, normalised weight of each trajectory
    ans <- list(log.likelihood=log.likelihood,traj=traj.particles,traj.weight=weight.particles/sum(weight.particles))

    return(ans)

}


# my_bootstrapParticleFilter <- function(fitmodel, n.particles)
# {

#     ############################################################################################
#     ## This function compute the marginal log.likelihood of the data (fitmodel$data) 
#     ## given the parameters (fitmodel$theta) using a particle filter
#     ############################################################################################

#     # Create some temporary variables (useful to avoid repetition of long names with $)
#     data <- fitmodel$data
#     theta <- fitmodel$theta

#     ############################################################################################
#     ## Initialisation of the algorithm
#     ############################################################################################

#     # Marginal log-likelihood is set to 0 and will be updated during the fitering steps
#     log.likelihood <- 0

#     # Initialise the particles:
#     # initial state: vector calculated by the function fitmodel$initialise.state evaluated at theta
#     initialise.particle  <- fitmodel$initialise.state(theta)

#     # particle states can be stored in a list
#     state.particles <- rep(list(initialise.particle),n.particles)
#     # weight: initially equal for all the particles 
#     # particle weight can be stored in a vector
#     weight.particles <- rep(1/n.particles,length=n.particles)

#     # initialise time variable
#     current.time <- 0

#     ############################################################################################
#     # Loop over observation time: resample, propagate, weight
#     ############################################################################################
#     for(i in seq_len(nrow(data))){

#         # extract next observtaion time from data
#         next.time <- data$time[i]

#         # resample particles according to their weight (normalization of the weight is done in the function sample())
#         index.resampled <- sample(x=n.particles,size=n.particles,replace=T,prob=weight.particles)
#         state.particles <- state.particles[index.resampled]

#         # propagate particles to the next observation time
#         for(p in 1:n.particles){
            
#             # extract current state of the particle 
#             current.state.particle <- unlist(state.particles[p])

#             # simulate from current observation time to next observation time
#             traj <- fitmodel$simulate.model(theta=theta,state.init=current.state.particle,times=c(current.time,next.time))

#             # compute particle weight
#             weight.particles[p] <- exp(fitmodel$log.likelihood( data=data, model.traj=traj, theta= theta))

#             # Update state of the p particle
#             # You need to take last row of the traj data frame which corresponds to the next observation time
#             # Also make sure to only take state variables (traj contains also a time variable)
#             state.particles[[p]] <- traj[2,fitmodel$state.variables]

#         }

#         # update time
#         current.time <- next.time

#         # Update marginal log-likelihood
#         log.likelihood <- log.likelihood + log(mean(weight.particles))

#     }


#     # Return marginal log-likelihood
#     return(log.likelihood)

# }
