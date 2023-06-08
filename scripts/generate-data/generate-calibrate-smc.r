library("here")
library("fitR")

data(models)

source(here::here("scripts", "solutions", "our-smc.r"))
source(here::here("scripts", "snippets", "calibrate-smc.r"))

dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)
save(calibrateSmc, file = file.path(dataDir, "calibrateSmc.rdata"))
