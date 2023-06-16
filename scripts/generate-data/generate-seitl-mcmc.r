library("fitR")

data(fluTdc1971)

dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)

setMcmcScript <- here::here("scripts", "snippets", "set-mcmc.r")
runMcmcScript <- here::here("scripts", "snippets", "run-mcmc.r")
thetaInitScript <- here::here("scripts", "snippets", "theta-init.r")
infoPriorScript <- here::here("scripts", "snippets", "seitl-info-prior.r")

source(setMcmcScript)
source(runMcmcScript)

mcmcSeitl <- my_mcmcTdc

source(thetaInitScript)
nIterations <- 50000

initTheta <- theta1
source(runMcmcScript)
mcmcSeitlTheta1 <- my_mcmcTdc

initTheta <- theta2
source(runMcmcScript)
mcmcSeitlTheta2 <- my_mcmcTdc

## informative priors
source(infoPriorScript)
seitlDeter$dPrior <- seitlInfoPrior

initTheta <- theta1
source(runMcmcScript)
mcmcSeitlInfoPriorTheta1 <- my_mcmcTdc

initTheta <- theta2
source(runMcmcScript)
mcmcSeitlInfoPriorTheta2 <- my_mcmcTdc

## SEIT4L model
tmp <- readLines(setMcmcScript)
tmp <- sub("my_fitmodel <- .*$", "my_fitmodel <- seit4lDeter", tmp)
tmp <- sub(
  "my_initState <- .*$",
  paste0(
    "my_initState <- ",
    "c(S = 279, E = 0, I = 2, T1 = 3, T2 = 0, T3 = 0, T4 = 0, L = 0, Inc = 0)"
  ), tmp)
source(textConnection(tmp))

nIterations <- 5000

source(runMcmcScript)
mcmcSeit4l <- my_mcmcTdc

source(thetaInitScript)
nIterations <- 50000

initTheta <- theta1
source(runMcmcScript)
mcmcSeit4lTheta1 <- my_mcmcTdc

initTheta <- theta2
source(runMcmcScript)
mcmcSeit4lTheta2 <- my_mcmcTdc

## informative priors
source(infoPriorScript)
seitlDeter$dPrior <- seitlInfoPrior

initTheta <- theta1
source(runMcmcScript)
mcmcSeit4lInfoPriorTheta1 <- my_mcmcTdc

initTheta <- theta2
source(runMcmcScript)
mcmcSeit4lInfoPriorTheta2 <- my_mcmcTdc

save(
  mcmcSeitl, mcmcSeit4l, file = here::here("data", "mcmcTdcDeterShortRun.rdata")
)

save(
  mcmcSeitlTheta1, mcmcSeitlTheta2,
  mcmcSeitlInfoPriorTheta1, mcmcSeitlInfoPriorTheta2,
  mcmcSeit4lTheta1, mcmcSeit4lTheta2,
  mcmcSeit4lInfoPriorTheta1, mcmcSeit4lInfoPriorTheta2,
  file = here::here("data", "mcmcTdcDeterLongRun.rdata")
)
