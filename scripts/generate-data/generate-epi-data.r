library("here")
library("fitR")

data(models)

parameters <- c(R_0 = 1.8, D_inf = 2)
initState <- c(S = 999, I = 1, R = 0)

epi1 <- rTrajObs(sirDeter, parameters, initState, 1:100)
epi1 <- epi1[1:38, ]
epi1 <- epi1[c("time", "obs")]

parameters <- c(R_0 = 2.5, D_inf = 2, RR = 0.1)
epi2 <- rTrajObs(sirReportingDeter, parameters, initState, 1:100)
epi2 <- epi1[1:26, ]
epi2 <- epi1[c("time", "obs")]

parameters <- c(R_0 = 2.1, D_inf = 3)
initState <- c(S = 999, I = 1, R = 0)

epi3 <- rTrajObs(sirDeter, parameters, initState, 1:100)
epi3 <- epi3[1:44, ]
epi3 <- epi3[c("time", "obs")]

parameters <- c(R_0 = 2.4, D_inf = 8, RR = 0.25)
initState <- c(S = 999, I = 1, R = 0)

epi4 <- rTrajObs(sirReportingDeter, parameters, initState, 1:100)
epi4 <- epi4[c("time", "obs")]

## save
dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)
save(epi1, epi2, epi3, epi4, file = file.path(dataDir, "epi.rdata"))
