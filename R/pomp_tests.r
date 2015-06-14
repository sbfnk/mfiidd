data(SEITL_pomp)
data(SEIT2L_pomp)

theta.bad.guess <- c(R0 = 2, D_lat = 2, D_inf = 2, alpha = 0.9, D_imm = 13, rho = 0.85)
init.state.bad.guess <- c(S = 250, E = 0, I = 4, T = 0, L = 30, Inc = 0)
init.state.bad.guess2 <- c(S = 250, E = 0, I = 4, T1 = 0, T2 = 0, L = 30, Inc = 0)

pomp.init.state.bad.guess <- init.state.bad.guess
pomp.init.state.bad.guess2 <- init.state.bad.guess2

names(pomp.init.state.bad.guess) <- paste(names(init.state.bad.guess), 0, sep = ".")
names(pomp.init.state.bad.guess2) <- paste(names(init.state.bad.guess2), 0, sep = ".")

## deterministic trajectory
traj <- trajectory(SEITL_pomp,
                   params = c(theta.bad.guess, pomp.init.state.bad.guess),
                   state = TRUE, obs = TRUE)
traj2 <- trajectory(SEIT2L_pomp,
                    params = c(theta.bad.guess, pomp.init.state.bad.guess),
                    state = TRUE, obs = TRUE)

## stochastic simulation
sim1 <- simulate(SEITL_pomp, params = c(theta.bad.guess, pomp.init.state.bad.guess))
sim2 <- simulate(SEIT2L_pomp, params = c(theta.bad.guess, pomp.init.state.bad.guess))

## trajectory matching
tm <- traj.match(SEITL_pomp,start = c(theta.bad.guess, pomp.init.state.bad.guess),
                 est = names(theta.bad.guess))
tm2 <- traj.match(SEIT2L_pomp,start=c(theta.bad.guess, pomp.init.state.bad.guess2),
                  est = names(theta.bad.guess))

tm_sim <- simulate(mf, nsim = 10, as.data.frame = TRUE, include.data = TRUE)
tm_sim2 <- simulate(mf2, nsim = 10, as.data.frame = TRUE, include.data = TRUE)

p <- ggplot(tm_sim, aes(x = time, y = obs, group = sim, alpha = (sim == "data")))
p <- p + scale_alpha_manual("", values = c(`TRUE` = 1, `FALSE` = 0.2),
                            labels = c(`FALSE` = "simulation",`TRUE` = "data"))
p <- p + geom_line()

p2 <- ggplot(tm2_sim, aes(x = time, y = obs, group = sim, alpha = (sim == "data")))
p2 <- p2 + scale_alpha_manual("", values = c(`TRUE` = 1, `FALSE` = 0.2),
                            labels = c(`FALSE` = "simulation",`TRUE` = "data"))
p2 <- p2 + geom_line()

## MIF
prop.sd <- theta.bad.guess
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
require(foreach)
require(doMC)
options(cores = 2)
registerDoMC()
mcopts <- list(set.seed = TRUE)
paropts <- list(.options.multicore = mcopts)

system.time(pm <- foreach (i=1:2, .combine=c, .options.multicore = mcopts) %dopar%
            { pmcmc(mf, Nmcmc = 3000, Np = 50, proposal = mvn.diag.rw(prop.sd),
                    verbose = TRUE, max.fail = Inf)
            })
traces <- conv.rec(pm,names(theta.bad.guess))
prop.sd <- apply(sapply(names(theta.bad.guess), function(x) sapply(traces[, x], sd)), 2, mean)
system.time(pm <- foreach (i=1:2, .combine=c, .options.multicore = mcopts) %dopar%
            { pmcmc(mf, Nmcmc = 3000, Np = 50, proposal = mvn.diag.rw(prop.sd),
                    verbose = TRUE, max.fail = Inf)
            })
plot(traces[, "R0"])
