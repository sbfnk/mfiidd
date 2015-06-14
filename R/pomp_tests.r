data(SEITL_pomp)
data(SEIT2L_pomp)

theta <-
    c(R0 = 2, D_lat = 2, D_inf = 2, alpha = 0.9, D_imm = 13, rho = 0.85)
SEITL.init.state <-
    c(S.0 = 250, E.0 = 0, I.0 = 4, T.0 = 0, L.0 = 30, Inc.0 = 0)
SEIT2L.init.state <-
    c(S.0 = 250, E.0 = 0, I.0 = 4, T1.0 = 0, T2.0 = 0, L.0 = 30, Inc.0 = 0)

## deterministic trajectory
traj <- trajectory(SEITL_pomp, params = c(theta, SEITL.init.state),
                   as.data.frame = TRUE)
traj2 <- trajectory(SEIT2L_pomp,
                    params = c(theta, SEIT2L.init.state),
                    state = TRUE, obs = TRUE)

## stochastic simulation
sim1 <- simulate(SEITL_pomp, params = c(theta, pomp.init.state))
sim2 <- simulate(SEIT2L_pomp, params = c(theta, pomp.init.state))

## trajectory matching
tm <- traj.match(SEITL_pomp,start = c(theta, pomp.init.state),
                 est = names(theta))
tm2 <- traj.match(SEIT2L_pomp,start=c(theta, pomp.init.state2),
                  est = names(theta))

tm_sim <- simulate(tm, nsim = 10, as.data.frame = TRUE, include.data = TRUE)
tm_sim2 <- simulate(tm2, nsim = 10, as.data.frame = TRUE, include.data = TRUE)

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

mf <- mif(tm, Nmif = 50, Np = 1000, cooling.fraction = 0.8, rw.sd = prop.sd)
mf2 <- mif(tm2, Nmif = 50, Np = 1000, cooling.fraction = 0.9, rw.sd = prop.sd)

mf_sim <- simulate(mf, nsim = 10, as.data.frame = TRUE, include.data = TRUE)
mf_sim2 <- simulate(mf2, nsim = 10, as.data.frame = TRUE, include.data = TRUE)

p <- ggplot(mf_sim, aes(x = time, y = obs, group = sim, alpha = (sim == "data")))
p <- p + scale_alpha_manual("", values = c(`TRUE` = 1, `FALSE` = 0.2),
                            labels = c(`FALSE` = "simulation",`TRUE` = "data"))
p <- p + geom_line()

p2 <- ggplot(mf_sim2, aes(x = time, y = obs, group = sim, alpha = (sim == "data")))
p2 <- p2 + scale_alpha_manual("", values = c(`TRUE` = 1, `FALSE` = 0.2),
                            labels = c(`FALSE` = "simulation",`TRUE` = "data"))
p2 <- p2 + geom_line()

## pMCMC
require('foreach')
require('doMC')
options(cores = 2)
registerDoMC()
mcopts <- list(set.seed = TRUE)
paropts <- list(.options.multicore = mcopts)

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

