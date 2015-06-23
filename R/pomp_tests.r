library('fitR')

example(SEITL_pomp)
example(SEIT4L_pomp)

theta <-
    c(R0 = 2, D_lat = 2, D_inf = 2, alpha = 0.9, D_imm = 13, rho = 0.85)
SEITL_init_state <-
    c(S.0 = 250, E.0 = 0, I.0 = 4, T.0 = 0, L.0 = 30, Inc.0 = 0)
SEIT4L_init_state <-
    c(S.0 = 250, E.0 = 0, I.0 = 4, T1.0 = 0, T2.0 = 0, L.0 = 30, Inc.0 = 0)

SEITL_tm <- traj.match(SEITL_pomp,start = c(theta, SEITL_init_state),
                       est = names(theta))
prop.sd <- theta
prop.sd[] <- 0.01
SEITL_pm <- pmcmc(SEITL_tm, Nmcmc = 100000, Np = 200, proposal = prop.sd, verbose = TRUE, max.fail = Inf, adapt.size.start = 100, adapt.shape.start = 5000, max.scaling = 50, adapt.size.cooling = 0.999)

SEITL_pm <- pmcmc(SEITL_tm, Nmcmc = 100000, Np = 200, proposal = prop.sd, verbose = TRUE, max.fail = Inf, adapt.size.start = 100, adapt.shape.start = 5000, max.scaling = 50, adapt.size.cooling = 0.999)

trace <- conv.rec(SEITL_pm, names(theta))
empirical.mode <- apply(trace, 2, median)
empirical.cov <- cov(trace)

SEITL_pm <- pmcmc(SEITL_pm, start = c(empirical.mode, SEITL_init_state), Nmcmc = 10000, Np = 200, proposal = empirical.cov, verbose = TRUE, max.fail = Inf, adapt.size.start = 100, max.scaling = 50, adapt.size.cooling = 0.999)


## deterministic trajectory
SEITL_traj <- trajectory(SEITL_pomp, params = c(theta, SEITL_init_state),
                         as.data.frame = TRUE)
SEITL4_traj2 <- trajectory(SEIT4L_pomp,
                           params = c(theta, SEIT4L_init_state),
                           state = TRUE, obs = TRUE)

## stochastic simulation
sim1 <- simulate(SEITL_pomp, params = c(theta, SEITL_init_state), states = TRUE)
sim2 <- simulate(SEIT4L_pomp, params = c(theta, SEIT4L_init_state))

## trajectory matching
SEITL_tm <- traj.match(SEITL_pomp,start = c(theta, SEITL_init_state),
                       est = names(theta))
SEIT4L_tm <- traj.match(SEIT4L_pomp,start=c(theta, SEIT4L_init_state),
                        est = names(theta))

SEITL_tm <- simulate(tm, nsim = 10, as.data.frame = TRUE, include.data = TRUE)
SEIT4L_tm <- simulate(tm2, nsim = 10, as.data.frame = TRUE, include.data = TRUE)

p <- ggplot(tm_sim, aes(x = time, y = obs, group = sim, alpha = (sim == "data")))
p <- p + scale_alpha_manual("", values = c(`TRUE` = 1, `FALSE` = 0.2),
                            labels = c(`FALSE` = "simulation",`TRUE` = "data"))
p <- p + geom_line()

p2 <- ggplot(tm2_sim, aes(x = time, y = obs, group = sim, alpha = (sim == "data")))
p2 <- p2 + scale_alpha_manual("", values = c(`TRUE` = 1, `FALSE` = 0.2),
                              labels = c(`FALSE` = "simulation",`TRUE` = "data"))
p2 <- p2 + geom_line()

## MIF
prop.sd <- theta
prop.sd[] <- 0.01

SEITL_mf <- mif(SEITL_tm, Nmif = 50, Np = 1000, cooling.fraction = 0.8, rw.sd = prop.sd)
SEIT4L_mf <- mif(SEIT4L_tm, Nmif = 50, Np = 1000, cooling.fraction = 0.9, rw.sd = prop.sd)

SEITL_mf_sim <- simulate(mf, nsim = 10, as.data.frame = TRUE, include.data = TRUE)
SEIT4L_mf_sim <- simulate(mf2, nsim = 10, as.data.frame = TRUE, include.data = TRUE)

p <- ggplot(mf_sim, aes(x = time, y = obs, group = sim, alpha = (sim == "data")))
p <- p + scale_alpha_manual("", values = c(`TRUE` = 1, `FALSE` = 0.2),
                            labels = c(`FALSE` = "simulation",`TRUE` = "data"))
p <- p + geom_line()

p2 <- ggplot(mf_sim2, aes(x = time, y = obs, group = sim, alpha = (sim == "data")))
p2 <- p2 + scale_alpha_manual("", values = c(`TRUE` = 1, `FALSE` = 0.2),
                              labels = c(`FALSE` = "simulation",`TRUE` = "data"))
p2 <- p2 + geom_line()

pf <- pfilter(SEITL_pomp, params = c(theta, SEITL_init_state), deterministic = TRUE)
pf <- pfilter(SEITL_pomp, params = c(theta, SEITL_init_state), Np = 100)


## pMCMC
require('foreach')
require('doMC')
options(cores = 2)
registerDoMC()
mcopts <- list(set.seed = TRUE)
paropts <- list(.options.multicore = mcopts)

prop.sd <- theta
prop.sd[] <- 0.01
SEITL_pm <- pmcmc(SEITL_tm, Nmcmc = 1000, Np = 200, proposal = prop.sd, verbose = TRUE, max.fail = Inf, adapt.size.start = 100, adapt.shape.start = 100)
trace <- conv.rec(SEITL_pm, names(theta))
## calculate the empirical covariance matrix
empirical.cov <- cov(log(trace))
## re-run pMCMC using the empirical correlations
SEITL_pm <- pmcmc(SEITL_tm, Nmcmc = 1000, Np = 200, proposal = mvn.rw(empirical.cov), 
                  verbose = TRUE, max.fail = Inf)
trace <- conv.rec(SEITL_pm, names(theta))
plot(trace)

require('coda')

system.time(pm50 <- foreach (i=1:2, .combine=c, .options.multicore = mcopts) %dopar%
            { pmcmc(mf, Nmcmc = 3000, Np = 50, proposal = mvn.diag.rw(prop.sd),
                    verbose = TRUE, max.fail = Inf)
            })
traces <- conv.rec(pm50,names(theta))
prop.sd.new <- apply(sapply(names(theta), function(x) sapply(traces[, x], sd)), 2, mean)

accRate <- 1 - rejectionRate(traces)

system.time(pm50 <- foreach (i=1:2, .combine=c, .options.multicore = mcopts) %dopar%
            { pmcmc(mf, Nmcmc = 3000, Np = 50, proposal = mvn.diag.rw(prop.sd.new),
                    verbose = TRUE, max.fail = Inf)
            })
traces <- conv.rec(pm50,names(theta))
rejectionRate(traces)
plot(traces[, "R0"])

system.time(pm400 <- foreach (i=1:2, .combine=c, .options.multicore = mcopts) %dopar%
            { pmcmc(mf, Nmcmc = 3000, Np = 400, proposal = mvn.diag.rw(prop.sd),
                    verbose = TRUE, max.fail = Inf)
            })
traces <- conv.rec(pm400,names(theta))
prop.sd.new <- apply(sapply(names(theta), function(x) sapply(traces[, x], sd)), 2, mean)

accRate <- 1 - rejectionRate(traces)

system.time(pm400 <- foreach (i=1:2, .combine=c, .options.multicore = mcopts) %dopar%
            { pmcmc(mf, Nmcmc = 3000, Np = 400, proposal = mvn.diag.rw(prop.sd.new),
                    verbose = TRUE, max.fail = Inf)
            })
traces <- conv.rec(pm400,names(theta))
rejectionRate(traces)
plot(traces)

pomp_to_fitmodel <- function(pompObject, stochastic = TRUE)
{
    if (stochastic)
    {
        ret <- function(theta, init.state, times)
        {
            pomp.init <- init.state
            names(pomp.init) <- paste(names(init.state), 0, sep = ".")
            sim <- simulate(pompObject, nsim = 1, params = c(theta, pomp.init), states = TRUE, as.data.frame = TRUE, times = times, t0 = times[1])
            sim$sim <- NULL
            return(sim)
        }
    } else
    {
        ret <- function(theta, init.state, times)
        {
            pomp.init <- init.state
            names(pomp.init) <- paste(names(init.state), 0, sep = ".")
            sim <- trajectory(pompObject, nsim = 1, params = c(theta, pomp.init), states = TRUE, as.data.frame = TRUE, times = times, t0 = times[1])
            sim$traj <- NULL
            return(sim)
        }
    }
    return(ret)
}

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
    state.particles  <- rep(list(init.state), n.particles)

    # Weight: initially equal for all the particles 
    # particle weight can be stored in a vector
    weight.particles <- rep(1/n.particles, length = n.particles)

    # Initialise time variable
    current.time <- 0

    ## Loop over observation times: resample, propagate, weight
    for(i in seq_len(nrow(data))){

        # Extract next data point (must be a vector)
        data.point <- unlist(data[i, ])
        next.time <- data.point["time"]

        # Resample particles according to their weights. 
        # You can use the `sample` function of R
        # (normalisation of the weights is done in the function)
        index.resampled <- sample(x = n.particles,
                                  size = n.particles,
                                  replace = TRUE,
                                  prob = weight.particles)
        state.particles <- state.particles[index.resampled]

        ## Loop over particles: propagate and weight
        for(p in 1:n.particles){

            # Extract current state of the particle 
            current.state.particle <- state.particles[[p]]

            # Propagate the particle from current observation time 
            # to the next one using the function `fitmodel$simulate`
            traj <- fitmodel$simulate(theta = theta,
                                      init.state = current.state.particle,
                                      times = c(current.time,next.time))

            # Extract state of the model at next observation time
            # Also make sure that model.point is a vector
            model.point <- unlist(traj[2,fitmodel$state.names])

            # Weight the particle with the likelihood of the observed 
            # data point using the function `fitmodel$dPointObs`
            weight.particles[p] <-
                exp(fitmodel$dPointObs(data.point = data.point,
                                       model.point = model.point,
                                       theta = theta))

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
