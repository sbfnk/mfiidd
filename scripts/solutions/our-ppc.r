# @knitr our_postPredCheck
# This is a function that takes 4 arguments:
# - trace, a data frame containing samples from the posterior
#   distribution, one column per parameter
# - nSamples, the number of samples to take
# - fitmodel, the model we use to generate replicates
# - initState, the initial state
# - data, the data set we have fit the model to
# It returns the two-sided p-value for the maximal observation
# in the data with respect to the model.
my_postPredCheck <- function(trace, nSamples, fitmodel, initState, data) {
  # calculate maximum in obs column of data
  maxData <- max(data$obs)

  # draw nSamples random numbers between 1
  # and nSamples using the `samples` function
  samples <- sample(seq_len(nrow(trace)), nSamples)

  # initialise vector of model maxima
  maxModel <- c()

  # loop over samples
  for (i in samples) {
    # get i'th column from the trace, unlist
    # (to convert to a vector) and assign to parameter
    # vector theta
    theta <- unlist(trace[i, ])

    # use rObsTraj to generate
    # observation trajectory using theta
    obsTraj <- rTrajObs(fitmodel, theta, initState, data$time)

    # calculate maximum in model and add to maxModel vector
    maxModel <- c(maxModel, max(obsTraj$obs))
  }

  # calculate quantiles of model maxima
  maxModelQuant <- quantile(maxModel, probs = c(0.025, 0.975))

  # calculate 2-sided p-value,
  # that is the proportion of elements of maxModel which are
  # either greater or equal or less or equal (whichever is
  # less) and  multiply by 2 (because it is a 2-sided test)
  pvalue <- min(
    sum(maxModel <= maxData),
    sum(maxModel >= maxData)
  ) / nSamples * 2

  # return two-sided p-value
  return(pvalue)
}
