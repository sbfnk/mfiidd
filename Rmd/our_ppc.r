# @knitr our_postPredCheck
# This is a function that takes 4 arguments:
# - trace, a data frame containing samples from the posterior
#   distribution, one column per parameter 
# - n.samples, the number of samples to take
# - fitmodel, the model we use to generate replicates
# - init.state, the initial state
# - data, the data set we have fit the model to
# It returns the two-sided p-value for the maximal observation
# in the data with respect to the model.
my_postPredCheck <- function(trace, n.samples, fitmodel, init.state, data) {

    # calculate maximum in obs column of data
    max.data <- max(data$obs)
    
    # draw n.samples random numbers between 1
    # and n.samples using the `samples` function 
    samples <- sample(seq_len(nrow(trace)), n.samples)

    # initialise vector of model maxima
    max.model <- c()
    
    # loop over samples
    for (i in samples) {

        # get i'th column from the trace, unlist
        # (to convert to a vector) and assign to parameter
        # vector theta
        theta <- unlist(trace[i, ])

        # use rObsTraj to generate
        # observation trajectory using theta
        obs.traj <- rTrajObs(fitmodel, theta, init.state, data$time)

        # calculate maximum in model and add to max.model vector
        max.model <- c(max.model, max(obs.traj$obs))
    }

    # calculate quantiles of model maxima
    max.model.quant <- quantile(max.model, probs = c(0.025, 0.975))

    # calculate 2-sided p-value,
    # that is the proportion of elements of max.model which are
    # either greater or equal or less or equal (whichever is
    # less) and  multiply by 2 (because it is a 2-sided test)
    pvalue <- min(sum(max.model <= max.data),
                  sum(max.model >= max.data)) / n.samples * 2
    
    # return two-sided p-value
    return(pvalue)
}
