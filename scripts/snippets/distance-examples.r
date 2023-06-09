# @knitr distance_examples
# mean absolute distance
ssMeanAbsDistance <- function(sumStats, dataObs, modelObs) {
  # calculate the absolute distance of each summary statistic,
  # and take the mean
  res <- mean(sapply(sumStats, function(x) {
    abs(x(modelObs) - x(dataObs))
  }))

  # return mean distance
  return(res)
}

# vector of absolute distances
ssAbsDistances <- function(sumStats, dataObs, modelObs) {
  # calculate the absolute distance of each summary statistic
  res <- sapply(sumStats, function(x) {
    abs(x(modelObs) - x(dataObs))
  })

  # set names of the vector of distances
  names(res) <- names(sumStats)

  # return vector of distance
  return(res)
}

# mean relative distance
ssMeanRelDistance <- function(sumStats, dataObs, modelObs) {
  # calculate the relative distance of each summary statistic,
  # and take the mean
  res <- mean(sapply(sumStats, function(x) {
    abs((x(modelObs) - x(dataObs)) / x(dataObs))
  }))

  # return mean distance
  return(res)
}

# vector of relative distances
ssRelDistances <- function(sumStats, dataObs, modelObs) {
  # calculate the relative distance of each summary statistic,
  res <- sapply(sumStats, function(x) {
    abs((x(obsTraj) - x(data)) / x(data))
  })

  # set names of the vector of distances
  names(res) <- names(sumStats)

  # return vector of distance
  return(res)
}
