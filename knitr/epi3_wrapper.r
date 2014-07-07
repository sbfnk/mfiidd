# @knitr epi3_wrapper
my_logPosterior_epi3 <- function(theta) {

    return(my_logPosterior(fitmodel = SIR,
                           theta = theta,
                           init.state = c(S = 999, I = 1, R = 0),
                           data = epi3))

}
