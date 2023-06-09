# @knitr epi3_wrapper
my_dLogPosteriorEpi3 <- function(theta) {
  return(my_dLogPosterior(
    fitmodel = SIR,
    theta = theta,
    initState = c(S = 999, I = 1, R = 0),
    data = epi3
  ))
}
