generate_html <- function(Rmd_name) {

	library(rmarkdown)
	library(knitr)

	repo_dir <- "/Users/Tonton/edu/Fit_course/mfiidd"
	input_dir <- file.path(repo_dir, "Rmd")
	output_dir <- file.path(repo_dir, "html")

	for(i in Rmd_name){
		render(file.path(input_dir,sprintf("%s.Rmd",i)),output_dir=output_dir,output_format=html_document(toc=TRUE, fig_width=5, fig_height=5, theme="flatly"))		
	}

}

main <- function() {

	# Rmd_name <- c("index","introduction","posterior_example","posterior_example_solution")
	# Rmd_name <- c("mcmc","mcmc_example","mcmc_example_solution","generate_samples")
	# Rmd_name <- c("mcmc_diagnostics","epi3_wrapper","mcmc_commands")
	# Rmd_name <- c("play_with_seitl","play_with_seitl_example")
	# Rmd_name <- c("mcmc_and_model_comparison","example_mcmc_SEITL","our_ppc","our_ppc_insert")
	# Rmd_name <- c("pmcmc","smc_example","smc_example_solution","pmcmc_solution")
	Rmd_name <- c("ABC","sumstat_examples","distance_examples","abc_solution")

	generate_html(Rmd_name)

}

main()