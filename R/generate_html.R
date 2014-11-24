generate_html <- function() {

	library(rmarkdown)
	library(knitr)

	render("/Users/Tonton/edu/Fit_course/mfiidd/Rmd_files/index.Rmd",output_format=html_document(toc=TRUE, fig_width=12, fig_height=8, theme="flatly"))

}

main <- function() {

	generate_html()

}

main()