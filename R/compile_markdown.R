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
# Rmd.files <- list.files(path.expand(paste0(mfiidd.dir, "/Rmd/")), ".*\\.Rmd$",
# 	full.names = FALSE)

# Rmd.files <- setdiff(Rmd.files, grep("pomp", Rmd.files, value=TRUE))

# for (file in Rmd.files) {

# 	Rmd.file <- path.expand(paste0(mfiidd.dir, "/", "Rmd/", file))
# 	html.file <- path.expand(paste0(mfiidd.dir, "/", "website/", sub("\\.Rmd$", ".html", file)))

# 	if (!file.exists(html.file) | file.mtime(Rmd.file) > file.mtime(html.file)) {
# 		render(Rmd.file, output_dir = path.expand(paste0(mfiidd.dir, "/website/")))
# 	}
# }


Rmd.files1 <- c("index","introduction","posterior_example","posterior_example_solution")
Rmd.files2 <- c("mcmc","mcmc_example","mcmc_example_solution","generate_samples")
Rmd.files3 <- c("mcmc_diagnostics","epi3_wrapper","mcmc_commands")
Rmd.files4 <- c("play_with_seitl","play_with_seitl_example")
Rmd.files5 <- c("mcmc_and_model_comparison","example_mcmc_SEITL","our_ppc","our_ppc_insert")
Rmd.files6 <- c("pmcmc","smc_example","smc_example_solution","pmcmc_solution")
Rmd.files7 <- c("pomp", "pomp_seitl_explanation")
Rmd.files <- c(Rmd.files7)

# Rmd.files <- c("ABC","sumstat_examples","distance_examples","abc_solution")


# Rmd.files <- path.expand(paste0(mfiidd.dir, "/Rmd/",c("","play_with_seitl_example"),".Rmd"))
# Rmd.files <- path.expand(paste0(mfiidd.dir, "/Rmd/",c("play_with_seitl","play_with_seitl_example"),".Rmd"))
# Rmd.files <- path.expand(paste0(mfiidd.dir, "/Rmd/",c("mcmc_and_model_comparison","example_mcmc_SEITL"),".Rmd"))

for(Rmd.file in Rmd.files){
	render(path.expand(paste0(mfiidd.dir, "/Rmd/", Rmd.file,".Rmd")), output_dir = path.expand(paste0(mfiidd.dir, "/website/"))) 
}

## Rmd.files <- path.expand(paste0(mfiidd.dir, "/Rmd/",c("play_with_seitl","play_with_seitl_example"),".Rmd"))

## for(Rmd.file in Rmd.files){
## 	render(Rmd.file, output_dir = path.expand(paste0(mfiidd.dir, "/website/"))) 
## }

## output_format=html_document(smart=FALSE, toc=TRUE, fig_width=5, fig_height=5, theme="spacelab", highlight="pygments")
