library(fitR)
example(SIR)
example(SIR_reporting)

parameters <- c(R0 = 1.8, D = 2)
state.init <- c(S = 999, I = 1, R = 0)

epi1 <- genObsTraj(SIR, parameters, state.init, 1:100)
epi1 <- epi1[1:38, ]
epi1 <- epi1[c("time", "obs")]

parameters <- c(R0 = 2.5, D = 2, RR = 0.1)
epi2 <- genObsTraj(SIR_reporting, parameters, state.init, 1:100)
epi2 <- epi1[1:26, ]
epi2 <- epi1[c("time", "obs")]

save(epi1, epi2, file="epi.rdata")

for (i in seq(1, 3, 0.1)) {
  parameters <- c(R0 = i, D = 2)
  post <- my_posterior(SIR, parameters, state.init, epi1)
  cat(i, post, "\n")
}
