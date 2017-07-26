start_me <- function(){

	library(rmarkdown)
	library(knitr)

	repo_dir <<- "/Users/Tonton/edu/Fit_course/mfiidd"
	rmd_dir <<- file.path(repo_dir, "Rmd")
	html_dir <<- file.path(repo_dir, "website")

	source(file.path(repo_dir, "R", "compile_mfiidd.R"))

}

create_data_for_ebola_project <- function() {

	dir_data <- file.path(html_dir, "data")
	
	if(!file.exists(dir_data)){
		dir.create(dir_data)
	}

	dbs_file <- "data_bases_synthesis.rds"
	dir_project <- "/Users/Tonton/work/projects/ebola"
	df_dbs <- readRDS(file.path(dir_project,"data","2014WestAfrica",dbs_file)) 

	# extract WHO data for L and SL

	df_data <- df_dbs %>% filter(report=="cases", type=="new", geo%in%c("Liberia","Sierra Leone","Guinea"), time=="weekly", source=="linelist_WHO", !is.na(value)) %>% 
	select(date, geo, value) %>% dplyr::rename(weekly_incidence=value)

	write.csv(df_data, file.path(dir_data,"Ebola_West_Africa_2014_WHO.csv"), row.names=FALSE)
	

	# extract data for Yambuku

	


}

generate_html <- function() {


	# Rmd_names <- grep(".Rmd",list.files(rmd_dir),value=TRUE)
	# Rmd_names <- c("index.Rmd")
	Rmd_names <- c("index.Rmd","ebola_project.Rmd")
	# Rmd_names <- c("index.Rmd","introduction.Rmd","posterior_example.Rmd","posterior_example_solution.Rmd")
	# Rmd_names <- c("mcmc.Rmd","mcmc_example.Rmd","mcmc_example_solution.Rmd","generate_samples.Rmd")
	# Rmd_names <- c("mcmc_diagnostics.Rmd","epi3_wrapper.Rmd","mcmc_commands.Rmd")
	# Rmd_names <- c("play_with_seitl.Rmd","play_with_seitl_example.Rmd")
	# Rmd_names <- c("mcmc_and_model_comparison.Rmd")#,"example_mcmc_SEITL.Rmd","our_ppc.Rmd","our_ppc_insert.Rmd")
	# Rmd_names <- c("pmcmc.Rmd","smc_example.Rmd","smc_example_solution.Rmd","pmcmc_solution.Rmd")
	# Rmd_names <- c("ABC.Rmd","sumstat_examples.Rmd","distance_examples.Rmd","abc_solution.Rmd")

	for(Rmd_name in Rmd_names){
		render(file.path(rmd_dir,Rmd_name),output_dir=html_dir,output_format=html_document(toc=TRUE, fig_width=5, fig_height=5, theme="flatly"))		
	}

}

main <- function() {

	start_me()
	# create_data_for_ebola_project()
	# generate_html()

	compile_mfiidd(practical = 4, clear.cache = FALSE)

}

main()