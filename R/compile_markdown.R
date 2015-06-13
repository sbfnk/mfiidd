mfiidd.dir <- switch(Sys.info()[["user"]],
	tonton = "~/edu/Fit_course/mfiidd",## Anton
	seb ="~/teaching/mfiidd"
	) 

library('rmarkdown')
library('coda')
library('plyr')
library('dplyr')
library('lattice')
library('ggplot2')
library('reshape2')

## first reinstall latest release of fitR so it is the one used in "render()"
install_github("sbfnk/fitR")

## recompile all html files
Rmd.files <- list.files(path.expand(paste0(mfiidd.dir, "/Rmd/")), ".*\\.Rmd$",
                        full.names = TRUE)

for (file in Rmd.files)
{
    render(file, output_dir = path.expand(paste0(mfiidd.dir, "/website/")), output_format=html_document(toc=TRUE, fig_width=5, fig_height=5, theme="flatly"))
}
