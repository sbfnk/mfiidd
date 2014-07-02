library(fitR)
example(SIR)
example(SIR_reporting)

parameters <- c(R0 = 1.8, D.inf = 2)
state.init <- c(S = 999, I = 1, R = 0)

epi1 <- genObsTraj(SIR, parameters, state.init, 1:100)
epi1 <- epi1[1:38, ]
epi1 <- epi1[c("time", "obs")]

parameters <- c(R0 = 2.5, D.inf = 2, RR = 0.1)
epi2 <- genObsTraj(SIR_reporting, parameters, state.init, 1:100)
epi2 <- epi1[1:26, ]
epi2 <- epi1[c("time", "obs")]

save(epi1, epi2, file="epi.rdata")

