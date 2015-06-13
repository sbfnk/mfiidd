mfiidd.dir <- "~/teaching/mfiidd" ## Seb

library('coda')
library('dplyr')
library('lattice')
library('ggplot2')

Rmd.files <- list.files(path.expand(paste0(mfiidd.dir, "/Rmd/")), ".*\\.Rmd$",
                        full.names = TRUE)

for (file in Rmd.files)
{
    render(file, output_dir = path.expand(paste0(mfiidd.dir, "/website/")))
}
