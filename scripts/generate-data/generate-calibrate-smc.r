library("here")
library("fitR")
library("future")

data(models)

my_particleFilter <- function(...) particleFilter(...)$margLogLike
source(here::here("scripts", "snippets", "calibrate-smc.r"))

dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)
save(calibrateSmc, file = file.path(dataDir, "calibrateSmc.rdata"))
