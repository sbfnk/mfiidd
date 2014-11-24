start_me <- function(){

	library(rmarkdown)
	library(knitr)

	repo_dir <<- "/Users/Tonton/edu/Fit_course/mfiidd"
	rmd_dir <<- file.path(repo_dir, "Rmd")
	html_dir <<- file.path(repo_dir, "website")


}

generate_html <- function() {


	Rmd_names <- grep(".Rmd",list.files(rmd_dir),value=TRUE)
	# Rmd_names <- c("index","introduction","posterior_example","posterior_example_solution")
	# Rmd_names <- c("mcmc","mcmc_example","mcmc_example_solution","generate_samples")
	# Rmd_names <- c("mcmc_diagnostics","epi3_wrapper","mcmc_commands")
	# Rmd_names <- c("play_with_seitl","play_with_seitl_example")
	# Rmd_names <- c("mcmc_and_model_comparison","example_mcmc_SEITL","our_ppc","our_ppc_insert")
	# Rmd_names <- c("pmcmc","smc_example","smc_example_solution","pmcmc_solution")
	# Rmd_names <- c("ABC","sumstat_examples","distance_examples","abc_solution")

	for(Rmd_name in Rmd_names){
		render(file.path(rmd_dir,Rmd_name),output_dir=html_dir,output_format=html_document(toc=TRUE, fig_width=5, fig_height=5, theme="flatly"))		
	}

}

main <- function() {

	start_me()
	generate_html()

}

main()