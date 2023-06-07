library("here")
library("purrr")
library("fitR")

modelDir <- here::here("scripts", "models")
dataFiles <- list.files(modelDir, pattern = "^example-")

safeSource <- purrr::safely(source)
## run twice to ensure dependencies are fulfilled
purrr::walk(dataFiles, \(x) safeSource(file.path(modelDir, x)))
purrr::walk(dataFiles, \(x) safeSource(file.path(modelDir, x)))

models <- gsub(
  "-", "_",
  sub("^example-(.+)\\.r", "\\1", dataFiles)
)

dataDir <- here::here("data")
dir.create(dataDir, showWarnings = FALSE)
save(list = models, file = here::here("data", "models.rdata"))
