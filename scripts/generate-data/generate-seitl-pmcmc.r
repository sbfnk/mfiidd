library("fitR")
library("furrr")
library("future")

data(fluTdc1971)
data(models)

dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)

my_particleFilter <- function(...) particleFilter(...)$margLogLike
setPmcmcScript <- here::here("scripts", "snippets", "set-pmcmc.r")
source(setPmcmcScript)

## informative priors
source(here::here("scripts", "snippets", "seitl-info-prior.r"))
seit4lStoch$dPrior <- seitlInfoPrior

cores <- future::availableCores()

future::plan(future::multicore)

nIterations <- 3000

start_time <- Sys.time()
source(here::here("scripts", "snippets", "run-parallel-pmcmc.r"))
end_time <- Sys.time()
pmcmcSeit4lInfoPrior128 <- my_Pmcmc
duration <- difftime(end_time, start_time, unit = "hours")
duration

## fewer particles
tmp <- readLines(setPmcmcScript)
tmp <- sub("my_nParticles <- .*$", "my_nParticles <- 16", tmp)
source(textConnection(tmp))

nIterations <- 3000

## re-run
start_time <- Sys.time()
source(here::here("scripts", "snippets", "run-parallel-pmcmc.r"))
end_time <- Sys.time()
duration <- difftime(end_time, start_time, unit = "hours")
pmcmcSeit4lInfoPrior16 <- my_Pmcmc
duration

save(
  pmcmcSeit4lInfoPrior16, pmcmcSeit4lInfoPrior128,
  file = here::here("data", "pmcmcSeit4lInfoPrior.rdata")
)
