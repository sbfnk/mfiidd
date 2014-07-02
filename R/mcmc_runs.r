my_logPosterior_epi1 <- function(theta) {

    return(my_logPosterior(fitmodel = SIR,
                        theta = theta,
                        state.init = c(S = 999, I = 1, R = 0),
                        data = epi1))

}

mcmc.run.R0 <- mcmcMH(target = my_logPosterior_epi1, theta.init = c(R0 = 3, D.inf = 2), proposal.sd = c(0.05, 0), n.iterations = 10000)

my_logPosterior_epi3 <- function(theta) {

    return(my_logPosterior(fitmodel = SIR,
                        theta = theta,
                        state.init = c(S = 999, I = 1, R = 0),
                        data = epi3))

}

mcmc.run.R0.Dinf <- mcmcMH(target = my_logPosterior_epi3, theta.init = c(R0 = 1, D.inf = 2), proposal.sd = c(0.01, 0.1), n.iterations = 10000)

save(mcmc.run.R0, mcmc.run.R0.Dinf, file = "mcmc.rdata")
