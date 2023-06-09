# @knitr our_postPredCheckInsert
# This is a function that takes 4 arguments:
# - trace, a data frame containing samples from the posterior
#   distribution, one column per parameter
# - nSamples, the number of samples to take
# - fitmodel, the model we use to generate replicates
# - initState, the initial state
# - data, the data set we have fit the model to
# It should return the two-sided p-value for the maximal observation
# in the data with respect to the model.
my_postPredCheck <- function(trace, nSamples, fitmodel, initState, data) {
  maxData <- # INSERT HERE: calculate maximum in obs column of data

    samples <- # INSERT HERE: draw nSamples random numbers between 1
    # and nSamples using the `samples` function

    # initialise vector of model maxima
    maxModel <- c()

  # loop over samples
  for (i in samples) {
    theta <- # INSERT HERE: get i'th column from the trace, unlist
      # (to convert to a vector) and assign to parameter
      # vector theta

      obsTraj <- # INSERT HERE: use rTrajObs to generate
      # observation trajectory using theta

      # calculate maximum in model and add to maxModel vector
      maxModel <- c(maxModel, max(obsTraj$obs))
  }

  pvalue <- # INSERT HERE: calculate 2-sided p-value,
    # that is the proportion of elements of maxModel which are
    # either greater or equal or less or equal (whichever is
    # less) and  multiply by 2 (because it is a 2-sided test)
    #

    # return two-sided p-value
    return(pvalue)
}
