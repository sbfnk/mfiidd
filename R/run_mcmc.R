library(fitR)

example(SEIT2L_stochastic)

if (length(args) > 1) {
    run <- args[2]
} else {
    run <- NULL
}

ssSize <- function(traj) {
    return(sum(traj$obs))
}

ssMax <- function(traj) {
    return(max(traj$obs))
}

ssMaxTime <- function(traj) {
    return(min(traj[which(traj$obs == max(traj$obs)), ]$time))
}

ssSum_13_24 <- function(traj) {
    return(sum(traj[traj$time > 12 & traj$time < 25, ]$obs))
}

ssMax_25_36 <- function(traj) {
    return(max(traj[traj$time > 24 & traj$time < 37, ]$obs))
}


ssSum_37_60 <- function(traj) {
    return(sum(traj[traj$time > 36, ]$obs))
}



ssAvgRelDistances <- function(sum.stats, fitmodel, theta, init.state, data) {

    obs.traj <- genObsTraj(fitmodel, theta, init.state, data$time)

    res <- mean(sapply(sum.stats, function(x) {
        abs((x(obs.traj) - x(data)) / x(data))
    }))

    return(res / length(sum.stats))

}

ABCLogPosterior <- function(epsilon, sum.stats, fitmodel, theta, init.state, data) {

    distance <- ssAvgRelDistances(sum.stats, fitmodel, theta, init.state, data)
    if (distance < epsilon) {
        log.density <- fitmodel$logPrior(theta)
    } else {
        log.density <- -Inf
    }

    return(list(log.density = log.density, distance = c(distance = distance)))
}

my_ABCLogPosterior_try_tdc <- function(theta) {

    ## init.state = c(S = 250, E = 0, I = 4, T1 = 0, T2 = 0, L = 30, Inc = 0)
        init.state <- c("S"=279,"E"=0,"I"=2,"T1"=3,"T2"=0,"L"=0,"Inc"=0)
    
    log.posterior <-
        ABCLogPosterior(0.02, list(ssSize = ssSize,
                                  ssMax = ssMax,
                                  ssMaxTime = ssMaxTime,
                                  ssSum_13_24 = ssSum_13_24,
                                  ssMax_25_36 = ssMax_25_36,
                                  ssSum_37_60 = ssSum_37_60),
                        SEIT2L_sto, theta, init.state, FluTdC1971)

    return(list(log.density = log.posterior$log.density,
                trace = c(theta, log.posterior$distance)))

}


init.theta.in4 <- c(R0 = 10.2095660, D.lat = 1.7725577, D.inf = 2.5733599, alpha = 0.6566896, D.imm = 19.3177172, rho = 0.8423282)

time.mcmc.abc.in4 <- system.time(mcmc.abc.in4 <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = init.theta.in4, n.iterations = 1000000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

if (is.null(run)) {
    saveRDS(time.mcmc.abc.in4, paste("mcmc_run.rds", sep = ""))
} else {
    saveRDS(time.mcmc.abc.in4, paste("mcmc_run_", run, ".rds", sep = ""))
}
