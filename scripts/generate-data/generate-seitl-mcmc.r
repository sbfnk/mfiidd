library("fitR")

data(fluTdc1971)

dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)

source(here::here("scripts", "snippets", "set-mcmc.r"))
source(here::here("scripts", "snippets", "run-mcmc.r"))

mcmcSeitl <- my_mcmcTdc
save(mcmcSeitl, file = here::here("data", "mcmcTdcDeterShortRun.rdata"))

source(here::here("scripts", "snippets", "theta-init.r"))
nIterations <- 100000

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

save(
  mcmcSeitlTheta1, mcmcSeitlTheta2,
  mcmcSeitlInfoPriorTheta1, mcmcSeitlInfoPriorTheta2,
  file = here::here("data", "mcmcTdcDeterLongRun.rdata")
)
