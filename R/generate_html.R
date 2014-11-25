start_me <- function(){

	library(rmarkdown)
	library(knitr)

	repo_dir <<- "/Users/Tonton/edu/Fit_course/mfiidd"
	rmd_dir <<- file.path(repo_dir, "Rmd")
	html_dir <<- file.path(repo_dir, "website")


}

generate_html <- function() {


	# Rmd_names <- grep(".Rmd",list.files(rmd_dir),value=TRUE)
	Rmd_names <- c("introduction.Rmd")
	# Rmd_names <- c("index.Rmd","introduction.Rmd","posterior_example.Rmd","posterior_example_solution.Rmd")
	# Rmd_names <- c("mcmc.Rmd","mcmc_example.Rmd","mcmc_example_solution.Rmd","generate_samples.Rmd")
	# Rmd_names <- c("mcmc_diagnostics.Rmd","epi3_wrapper.Rmd","mcmc_commands.Rmd")
	# Rmd_names <- c("play_with_seitl.Rmd","play_with_seitl_example.Rmd")
	# Rmd_names <- c("mcmc_and_model_comparison.Rmd","example_mcmc_SEITL.Rmd","our_ppc.Rmd","our_ppc_insert.Rmd")
	# Rmd_names <- c("pmcmc.Rmd","smc_example.Rmd","smc_example_solution.Rmd","pmcmc_solution.Rmd")
	# Rmd_names <- c("ABC.Rmd","sumstat_examples.Rmd","distance_examples.Rmd","abc_solution.Rmd")

	for(Rmd_name in Rmd_names){
		render(file.path(rmd_dir,Rmd_name),output_dir=html_dir,output_format=html_document(toc=TRUE, fig_width=5, fig_height=5, theme="flatly"))		
	}

}

main <- function() {

	start_me()
	generate_html()

}

main()