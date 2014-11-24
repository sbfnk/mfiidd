# @knitr distance_examples
# mean absolute distance
ssMeanAbsDistance <- function(sum.stats, data.obs, model.obs) {

    # calculate the absolute distance of each summary statistic,
    # and take the mean
    res <- mean(sapply(sum.stats, function(x) {
        abs(x(model.obs) - x(data.obs))
    }))

    # return mean distance
    return(res)

}

# vector of absolute distances
ssAbsDistances <- function(sum.stats, data.obs, model.obs) {

    # calculate the absolute distance of each summary statistic
    res <- sapply(sum.stats, function(x) {
        abs(x(model.obs) - x(data.obs))
    })

    # set names of the vector of distances
    names(res) <- names(sum.stats)

    # return vector of distance
    return(res)

}

# mean relative distance
ssMeanRelDistance <- function(sum.stats, data.obs, model.obs) {

    # calculate the relative distance of each summary statistic,
    # and take the mean 
    res <- mean(sapply(sum.stats, function(x) {
        abs((x(model.obs) - x(data.obs)) / x(data.obs))
    }))

    # return mean distance
    return(res)

}

# vector of relative distances
ssRelDistances <- function(sum.stats, data.obs, model.obs) {

    # calculate the relative distance of each summary statistic,
    res <- sapply(sum.stats, function(x) {
        abs((x(obs.traj) - x(data)) / x(data))
    })

    # set names of the vector of distances
    names(res) <- names(sum.stats)

    # return vector of distance
    return(res)

}

