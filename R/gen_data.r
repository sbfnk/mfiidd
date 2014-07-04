library(fitR)
example(SIR)
example(SIR_reporting)

parameters <- c(R0 = 1.8, D.inf = 2)
init.state <- c(S = 999, I = 1, R = 0)

epi1 <- genObsTraj(SIR, parameters, init.state, 1:100)
epi1 <- epi1[1:38, ]
epi1 <- epi1[c("time", "obs")]

parameters <- c(R0 = 2.5, D.inf = 2, RR = 0.1)
epi2 <- genObsTraj(SIR_reporting, parameters, init.state, 1:100)
epi2 <- epi1[1:26, ]
epi2 <- epi1[c("time", "obs")]

parameters <- c(R0 = 2.1, D.inf = 3)
init.state <- c(S = 999, I = 1, R = 0)

epi3 <- genObsTraj(SIR, parameters, init.state, 1:100)
epi3 <- epi3[1:44, ]
epi3 <- epi3[c("time", "obs")]

parameters <- c(R0 = 2.4, D.inf = 8, RR = 0.25)
init.state <- c(S = 999, I = 1, R = 0)

epi4 <- genObsTraj(SIR_reporting, parameters, init.state, 1:100)
epi4 <- epi4[c("time", "obs")]

save(epi1, epi2, epi3, epi4, file="epi.rdata")

