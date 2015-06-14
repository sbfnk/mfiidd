mfiidd.dir <- switch(Sys.info()[["user"]],
	Tonton = "~/edu/Fit_course/mfiidd",## Anton
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
# install_github("sbfnk/fitR")

## recompile all html files
Rmd.files <- list.files(path.expand(paste0(mfiidd.dir, "/Rmd/")), ".*\\.Rmd$",
	full.names = FALSE)


for (file in Rmd.files) {
	Rmd.file <- path.expand(paste0(mfiidd.dir, "/", "Rmd/", file))
	html.file <- path.expand(paste0(mfiidd.dir, "/", "website/",
		sub("\\.Rmd$", ".html", file)))
	if (!file.exists(html.file) | file.mtime(Rmd.file) > file.mtime(html.file)) {
		render(file, output_dir = path.expand(paste0(mfiidd.dir, "/website/")), output_format=html_document(smart=FALSE, toc=TRUE, fig_width=5, fig_height=5, theme="spacelab", highlight="pygments"))
	}
}
