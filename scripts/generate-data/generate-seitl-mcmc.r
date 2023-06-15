library("fitR")

data(fluTdc1971)

dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)

source(here::here("scripts", "snippets", "set-mcmc.r"))
source(here::here("scripts", "snippets", "run-mcmc.r"))

mcmcSeitl <- my_mcmcTdc

source(here::here("scripts", "snippets", "theta-init.r"))
nIterations <- 50000

initTheta <- theta1
source(here::here("scripts", "snippets", "run-mcmc.r"))
mcmcSeitlTheta1 <- my_mcmcTdc

initTheta <- theta2
source(here::here("scripts", "snippets", "run-mcmc.r"))
mcmcSeitlTheta2 <- my_mcmcTdc

## informative priors
source(here::here("scripts", "snippets", "seitl-info-prior.r"))
seitlDeter$dPrior <- seitlInfoPrior

initTheta <- theta1
source(here::here("scripts", "snippets", "run-mcmc.r"))
mcmcSeitlInfoPriorTheta1 <- my_mcmcTdc

initTheta <- theta2
source(here::here("scripts", "snippets", "run-mcmc.r"))
mcmcSeitlInfoPriorTheta2 <- my_mcmcTdc

## SEIT4L model
tmp <- readLines(setPmcmcScript)
tmp <- sub("my_fitmodel <- .*$", "my_fitmodel <- seit4lDeter", tmp)
tmp <- sub(
  "my_initState <- .*$",
  paste0(
    "my_initState <- ",
    "c(S = 279, E = 0, I = 2, T1 = 3, T2 = 0, T3 = 0, T4 = 0, L = 0, Inc = 0)"
  ), tmp)
source(textConnection(tmp))

source(here::here("scripts", "snippets", "set-mcmc.r"))
source(here::here("scripts", "snippets", "run-mcmc.r"))

mcmcSeit4l <- my_mcmcTdc

source(here::here("scripts", "snippets", "theta-init.r"))
nIterations <- 50000

initTheta <- theta1
source(here::here("scripts", "snippets", "run-mcmc.r"))
mcmcSeit4lTheta1 <- my_mcmcTdc

initTheta <- theta2
source(here::here("scripts", "snippets", "run-mcmc.r"))
mcmcSeit4lTheta2 <- my_mcmcTdc

## informative priors
source(here::here("scripts", "snippets", "seitl-info-prior.r"))
seitlDeter$dPrior <- seitlInfoPrior

initTheta <- theta1
source(here::here("scripts", "snippets", "run-mcmc.r"))
mcmcSeit4lInfoPriorTheta1 <- my_mcmcTdc

initTheta <- theta2
source(here::here("scripts", "snippets", "run-mcmc.r"))
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
