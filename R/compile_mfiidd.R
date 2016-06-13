##' Compiles mfiidd practicals
##'
##' @param Rmd.files character vector. Rmd files to compile (with or without the ".Rmd" extension). If not given, will compile all files ending on \code{.Rmd} in the \code{Rmd} directory
##' @param exclude character vector. Any files to exclude (with or without the ".Rmd" extension). This overrides any information given by \code{Rmd.files}
##' @param mfiidd.dir character vector. The directory in which \code{mfiidd} is located. Uses a default for Anton and Seb.
##' @param reinstall.fitR boolean. Whether fitR should be reinstalled (from github)
##' @param practical integer vector. If given, will only compile the files for the given practical session. This respects any information given in \code{Rmd.files} or \code{exclude}
##' @param clear.cache boolean. Whether to clear the cache before compiling
##' @author seb
compile_mfiidd <- function(Rmd.files = NULL, exclude = NULL, mfiidd.dir = NULL, reinstall.fitR = FALSE, practical = NULL, clear.cache = FALSE)
{

    require('rmarkdown')
    require('coda')
    require('plyr')
    require('dplyr')
    require('lattice')
    require('ggplot2')
    require('reshape2')

    if (is.null(mfiidd.dir))
    {
        mfiidd.dir <- switch(Sys.info()[["user"]],
           Tonton = "~/edu/Fit_course/mfiidd", ## Anton
           seb ="~/teaching/mfiidd" ## Seb
           )
    }

    if (reinstall.fitR)
    {
        ## first reinstall latest release of fitR so it is the one used in "render()"
        install_github("sbfnk/fitR")
    }

    if (is.null(Rmd.files))
    {
        Rmd.files <- list.files(path.expand(paste0(mfiidd.dir, "/Rmd/")), ".*\\.Rmd$",
            full.names = FALSE)
    }

    Rmd.files <- sub(".Rmd$", "", Rmd.files)

    if (!is.null(practical))
    {
        practical.files <- list(c("index","introduction","posterior_example","posterior_example_solution"), 
            c("mcmc","mcmc_example","mcmc_example_solution","generate_samples"), 
            c("mcmc_diagnostics","epi3_wrapper","mcmc_commands"),
            c("play_with_seitl","play_with_seitl_example"),
            c("mcmc_and_model_comparison","example_mcmc_SEITL","our_ppc","our_ppc_insert"), 
            c("pmcmc","smc_example","smc_example_solution","pmcmc_solution"), 
            c("pomp", "pomp_seitl_explanation"), 
            c("ABC","sumstat_examples","distance_examples","abc_solution"))

        practical <- as.integer(practical)

        if (any(is.na(practical)))
        {
            stop(sQuote("practical"), " must be a vector of interger")
        } else
        {
            Rmd.files <- unlist(practical.files[practical])
        }
    }

    if (!is.null(exclude))
    {
        exclude <- sub(".Rmd$", "", Rmd.files)
        Rmd.files <- setdiff(Rmd.files, exclude)
    }

    if (clear.cache)
    {
        unlink(paste0(mfiidd.dir, "/Rmd/cache"), recursive = TRUE)
    }

    if (length(Rmd.files) > 0) {

        Rmd_with_toc_float <- c("introduction", "mcmc", "mcmc_diagnostics", "play_with_seitl", "play_with_seitl_example", "mcmc_and_model_comparison", "example_mcmc_SEITL", "pmcmc", "pmcmc_solution")

        for (Rmd.file in Rmd.files) {
            render(
                path.expand(paste0(mfiidd.dir, "/Rmd/", Rmd.file,".Rmd")), 
                output_dir = path.expand(paste0(mfiidd.dir, "/website/")), 
                output_format = html_document(
                    toc = Rmd.file %in% Rmd_with_toc_float, 
                    toc_float = Rmd.file %in% Rmd_with_toc_float,
                    number_sections = TRUE,
                    theme = "yeti",
                    highlight = "tango",
                    smart = FALSE,
                    fig_width = 5,
                    fig_height = 5
                    )
                )
        }
    } else
    {
        stop("No files to compile")
    }
}

