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

	Rmd_name <- c("index","introduction")

	generate_html(Rmd_name)

}

main()