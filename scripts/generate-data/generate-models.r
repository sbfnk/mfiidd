library("here")
library("purrr")
library("fitR")

modelDir <- here::here("scripts", "models")
dataFiles <- list.files(modelDir, pattern = "^example-")

safeSource <- purrr::safely(source)
## run twice to ensure dependencies are fulfilled
purrr::walk(dataFiles, \(x) safeSource(file.path(modelDir, x)))
purrr::walk(dataFiles, \(x) source(file.path(modelDir, x)))

## convert to camel
models <- sub("^example-(.+)\\.r", "\\1", dataFiles)
models <- sub("^([A-Z0-9]+)\\-", "\\L\\1-", models, perl = TRUE)
models <- gsub("-([a-z0-9])", "\\U\\1", models, perl = TRUE)

dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)
save(list = models, file = here::here("data", "models.rdata"))
