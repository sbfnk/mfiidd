my_logPosterior_epi1 <- function(theta) {

    return(my_logPosterior(fitmodel = SIR,
                        theta = theta,
                        init.state = c(S = 999, I = 1, R = 0),
                        data = epi1))

}

## mcmc.run.R0 <- mcmcMH(target = my_logPosterior_epi1, init.theta = c(R0 = 3, D.inf = 2), proposal.sd = c(0.05, 0), n.iterations = 10000)
mcmc.epi1 <- mcmcMH(target = my_logPosterior_epi1, init.theta = c(R0 = 3, D.inf = 2), proposal.sd = c(0.05, 0), n.iterations = 10000)

my_logPosterior_epi3 <- function(theta) {

    return(my_logPosterior(fitmodel = SIR,
                        theta = theta,
                        init.state = c(S = 999, I = 1, R = 0),
                        data = epi3))

}

## mcmc.run.R0.Dinf <- mcmcMH(target = my_logPosterior_epi3, init.theta = c(R0 = 1, D.inf = 2), proposal.sd = c(0.01, 0.1), n.iterations = 10000)
mcmc.epi3 <- mcmcMH(target = my_logPosterior_epi3, init.theta = c(R0 = 1, D.inf = 2), proposal.sd = c(0.01, 0.1), n.iterations = 1000)

my_logPosterior_epi4 <- function(theta) {

    return(my_logPosterior(fitmodel = SIR_reporting,
                        theta = theta,
                        init.state = c(S = 999, I = 1, R = 0),
                        data = epi4))

}

## mcmc.run.R0.Dinf.RR <- mcmcMH(target = my_logPosterior_epi4, init.theta = c(R0 = 1, D.inf = 2, RR = 1), proposal.sd = c(0.01, 0.1, 0.01), n.iterations = 10000)
mcmc.epi4 <- mcmcMH(target = my_logPosterior_epi4, init.theta = c(R0 = 1, D.inf = 2, RR = 1), proposal.sd = c(0.01, 0.1, 0.01), n.iterations = 10000)

save(mcmc.epi1, mcmc.epi3, mcmc.epi4, file = "mcmc.rdata")

SEIT2L_sto$logPrior <- function(theta) {

	log.prior.R0 <- dunif(theta[["R0"]], min = 1, max = 50, log = TRUE)
	log.prior.latent.period <- dnorm(theta[["D.lat"]], mean = 2, log = TRUE)
	log.prior.infectious.period <- dnorm(theta[["D.inf"]], mean = 2, log = TRUE)
	log.prior.temporary.immune.period <- dunif(theta[["D.imm"]], min = 0, max = 50, log = TRUE)
	log.prior.probability.long.term.immunity <- dunif(theta[["alpha"]], min = 0, max = 1, log = TRUE)
	log.prior.reporting.rate <- dunif(theta[["rho"]], min = 0, max = 1, log = TRUE)
	
	return(log.prior.R0 + log.prior.latent.period + log.prior.infectious.period + log.prior.temporary.immune.period + log.prior.probability.long.term.immunity + log.prior.reporting.rate)

}

time.mcmc.abc.trial <- system.time(mcmc.abc.trial <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = theta, n.iterations = 10000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1)), proposal.sd = rep(1, 6)))

epsilon <- unname(quantile(mcmc.abc.trial$trace$distance, probs = 0.01))
init.theta.further <- unlist(mcmc.abc.trial$trace[which.min(mcmc.abc.trial$trace$distance), 1:6])

time.mcmc.abc.further <- system.time(mcmc.abc.further <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = init.theta.further, proposal.sd = init.theta.further/10, n.iterations = 10000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

theta.abc <- summary(mcmc(mcmc.abc.further$trace))$statistics[, 1]

epsilon.further <- unname(quantile(mcmc.abc.further$trace$distance, probs = 0.01))
init.theta.further2 <- unlist(mcmc.abc.further$trace[which.min(mcmc.abc.further$trace$distance), 1:6])

time.mcmc.abc.further2 <- system.time(mcmc.abc.further2 <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = init.theta.further, proposal.sd = init.theta.further/10, n.iterations = 100000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

epsilon.further2 <- unname(quantile(mcmc.abc.further2$trace$distance, probs = 0.1))
init.theta.further3 <- unlist(mcmc.abc.further2$trace[which.min(mcmc.abc.further2$trace$distance), 1:6])

time.mcmc.abc.further3 <- system.time(mcmc.abc.further3 <- mcmcMH(target = my_ABClogPosterior_try_tdc, init.theta = init.theta.further3, proposal.sd = init.theta.further3/10, n.iterations = 100000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

time.mcmc.abcvec.trial <- system.time(mcmc.abcvec.trial <- mcmcMH(target = my_ABCVecLogPosterior_try_tdc, init.theta = theta, n.iterations = 10000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1)), proposal.sd = rep(1, 6)))

epsilon <- unname(quantile(mcmc.abcvec.trial$trace$distance, probs = 0.01))
init.theta.further <- unlist(mcmc.abcvec.trial$trace[which.min(mcmc.abcvec.trial$trace$distance), 1:6])

time.mcmc.abcvec.further <- system.time(mcmc.abcvec.further <- mcmcMH(target = my_ABClogPosterior_try_tdc, init.theta = init.theta.further, proposal.sd = init.theta.further/10, n.iterations = 10000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

epsilon.further <- unname(quantile(mcmc.abcvec.further$trace$distance, probs = 0.01))
init.theta.further2 <- unlist(mcmc.abcvec.further$trace[which.min(mcmc.abcvec.further$trace$distance), 1:6])

time.mcmc.abcvec.further2 <- system.time(mcmc.abcvec.further2 <- mcmcMH(target = my_ABClogPosterior_try_tdc, init.theta = init.theta.further2, proposal.sd = init.theta.further2/10, n.iterations = 10000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

epsilon.further2 <- unname(quantile(mcmc.abcvec.further2$trace$distance, probs = 0.1))
init.theta.further3 <- unlist(mcmc.abcvec.further2$trace[which.min(mcmc.abcvec.further2$trace$distance), 1:6])

time.mcmc.abcvec.further3 <- system.time(mcmc.abcvec.further3 <- mcmcMH(target = my_ABClogPosterior_try_tdc, init.theta = init.theta.further3, proposal.sd = init.theta.further3/10, n.iterations = 100000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

time.mcmc.abc.in.trial <- system.time(mcmc.abc.in.trial <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = theta, n.iterations = 10000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1)), proposal.sd = rep(1, 6)))

epsilon.in <- unname(quantile(mcmc.abc.in.trial$trace$distance, probs = 0.01))
init.theta.in.further <- unlist(mcmc.abc.in.trial$trace[which.min(mcmc.abc.in.trial$trace$distance), 1:6])

time.mcmc.abc.in.further <- system.time(mcmc.abc.in.further <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = init.theta.in.further, proposal.sd = init.theta.in.further/10, n.iterations = 10000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

epsilon.in.further <- unname(quantile(mcmc.abc.in.further$trace$distance, probs = 0.01))
init.theta.in2 <- unlist(mcmc.abc.in.further$trace[which.min(mcmc.abc.in.further$trace$distance), 1:6])

time.mcmc.abc.in2 <- system.time(mcmc.abc.in2 <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = init.theta.in.further, proposal.sd = init.theta.in.further/10, n.iterations = 10000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

epsilon.in2 <- unname(quantile(mcmc.abc.in2$trace$distance, probs = 0.01))
init.theta.in3 <- unlist(mcmc.abc.in.further$trace[which.min(mcmc.abc.in.further$trace$distance), 1:6])

time.mcmc.abc.in3 <- system.time(mcmc.abc.in3 <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = init.theta.in3, n.iterations = 100000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

epsilon.in3 <- unname(quantile(mcmc.abc.in3$trace$distance, probs = 0.1))
init.theta.in4 <- unlist(mcmc.abc.in3$trace[which.min(mcmc.abc.in3$trace$distance), 1:6])

time.mcmc.abc.in4 <- system.time(mcmc.abc.in4 <- mcmcMH(target = my_ABCLogPosterior_try_tdc, init.theta = init.theta.in4, n.iterations = 1000000, limits = list(lower = c(R0 = 1, D.lat = 0, D.inf = 0, D.imm = 0, alpha = 0, rho = 0), upper = c(R0 = Inf, D.lat = Inf, D.inf = Inf, D.imm = Inf, alpha = 1, rho = 1))))

SEITL_pomp
