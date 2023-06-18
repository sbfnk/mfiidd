seitlSkelR <- function(x, t, params, ...) {

  with(as.list(c(x, params)), {
    beta <- R0 / D_inf
    epsilon <- 1 / D_lat
    nu <- 1 / D_inf
    tau <- 1 / D_imm
    N <- S + E + I + T + L

    trans <- c(
      beta * I / N * S,
      epsilon * E,
      nu * I,
      alpha * tau * T,
      (1 - alpha) * tau * T
    )

    DS <- -trans[1] + trans[5]
    DE <- trans[1] - trans[2]
    DI <- trans[2] - trans[3]
    DT <- trans[3] - trans[4] - trans[5]
    DL <- trans[4]
    DInc <- trans[2]

    c(DS, DE, DI, DT, DL, DInc)
  })
}

seitlSimR <- function(x, t, params, delta.t, ...) {
  with(as.list(c(x, params)), {
    beta = R0 / D_inf
    epsilon = 1 / D_lat
    nu = 1 / D_inf
    tau = 1 / D_imm

    N = S + E + I + T + L

    rate <- c(
      beta * I / N,
      epsilon,
      nu,
      alpha * tau,
      (1 - alpha) * tau
    )

    dN <- c(
      reulermultinom(n = 1, size = S, rate = rate[1], dt = delta.t),
      reulermultinom(n = 1, size = E, rate = rate[2], dt = delta.t),
      reulermultinom(n = 1, size = I, rate = rate[3], dt = delta.t),
      reulermultinom(n = 2, size = T, rate = rate[4:5], dt = delta.t)
    )

    S <- S - dN[1] + dN[5]
    E <- E + dN[1] - dN[2]
    I <- I + dN[2] - dN[3]
    T <- T + dN[3] - dN[4] - dN[5]
    L <- L + dN[4]
    Inc <- Inc + dN[2]

    c(S, E, I, T, L, Inc)
  })
}

seitlGenObsPointR <- function(x, t, params, ...) {
  with(as.list(c(x, params)), {
    c(obs = rpois(n = 1, ifelse(rho * Inc > 0, rho * Inc, 0)))
  })
}

seitlPointLikeR <- function(y, x, t, params, log, ...) {
  with(as.list(c(y, x, params)), {
    dpois(obs, ifelse(rho * Inc > 0, rho * Inc, 0), log)
  })
}

seitlPriorR <- function(params, log = FALSE, ...) {
  with(as.list(c(params)), {
    lik <- dunif(R0, 1, 50, log = TRUE) +
      dunif(D_lat, 0, 10, log = TRUE) +
      dunif(D_inf, 0, 15, log = TRUE) +
      dunif(D_imm, 0, 50, log = TRUE) +
      dunif(alpha, 0, 1, log = TRUE) +
      dunif(rho, 0, 1, log = TRUE)

    ifelse(log, lik, exp(lik))
  })
}

## construct pomp object
seitlPompR <- pomp(
  data = fluTdc1971[, c("time", "obs")],
  skeleton = vectorfield(seitlSkelR),
  rprocess = euler(step.fun = seitlSimR, delta.t = 0.1),
  rmeasure = seitlGenObsPointR,
  dmeasure = seitlPointLikeR,
  dprior = seitlPriorR,
  partrans = parameter_trans(
    log = c("R0", "D_inf", "D_lat", "D_imm", "alpha", "rho")
  ),
  times = "time",
  t0 = 1,
  accumvars = "Inc",
  paramnames = c("R0", "D_inf", "D_lat", "D_imm", "alpha", "rho"),
  statenames = c("S", "E", "I", "T", "L", "Inc"),
  obsnames = c("obs")
)
