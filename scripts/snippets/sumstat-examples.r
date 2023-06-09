# @knitr sumstat_examples
# maximum incidence before day 13
ssMax_0_12 <- function(traj) {
  return(max(traj[traj$time < 13, ]$obs, na.rm = T))
}

# maximum incidence between day 13 and day 24
ssMax_13_24 <- function(traj) {
  return(max(traj[traj$time > 12 & traj$time < 25, ]$obs, na.rm = T))
}

# maximum incidence between day 25 and day 36
ssMax_25_36 <- function(traj) {
  return(max(traj[traj$time > 24 & traj$time < 37, ]$obs, na.rm = T))
}

# maximum incidence after day 36
ssMax_37_60 <- function(traj) {
  return(max(traj[traj$time > 36, ]$obs, na.rm = T))
}

# cumulative incidence before day 13
ssSum_0_12 <- function(traj) {
  return(sum(traj[traj$time < 13, ]$obs, na.rm = T))
}

# cumulative incidence between day 13 and day 24
ssSum_13_24 <- function(traj) {
  return(sum(traj[traj$time > 12 & traj$time < 25, ]$obs, na.rm = T))
}

# cumulative incidence between day 25 and day 36
ssSum_25_36 <- function(traj) {
  return(sum(traj[traj$time > 24 & traj$time < 37, ]$obs, na.rm = T))
}

# cumulative incidence after day 36
ssSum_37_60 <- function(traj) {
  return(sum(traj[traj$time > 36, ]$obs, na.rm = T))
}

# maximum incidence along the whole trajectory
ssMax <- function(traj) {
  return(max(traj$obs, na.rm = T))
}

# timing of the epidemic peak
ssMaxTime <- function(traj) {
  return(min(traj[which(traj$obs == max(traj$obs)), ]$time, na.rm = T))
}

# final size of the epidemic
ssSize <- function(traj) {
  return(sum(traj$obs, na.rm = T))
}
