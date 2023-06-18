model seitlDeter {
  const k_erlang = 1
  const N = 1000
  const timestep = 1
  dim k(k_erlang)
  state S, E, I, T[k], L, Inc
  param R_0, D_lat, D_inf, alpha, D_imm, rho
  noise infection
  noise incubation
  noise loss_infectiousness
  noise immunity[k]
  noise loss_immunity
  obs Cases
  sub initial {
    E <- 0
    I <- 2
    T[k] <- (k == 0 ? 3 : 0)
    L <- 0
    Inc <- 0
    S <- N - I - T[0]
  }
  sub parameter {
    R_0 ~ uniform(1, 50)
    D_lat ~ uniform(0, 10)
    D_inf ~ uniform(0, 15)
    D_imm ~ uniform(0, 50)
    alpha ~ uniform(0, 1)
    rho ~ uniform(0, 1)
  }
  sub transition (delta=timestep) {
    inline beta = R_0/D_inf
    inline epsilon = 1/D_lat
    inline nu = 1/D_inf
    inline tau = 1/D_imm
    Inc <- 0
    infection ~ binomial(S, 1 - exp(-beta * I/N * timestep))
    incubation ~ binomial(E, 1 - exp(-epsilon * timestep))
    loss_infectiousness ~ binomial(I, 1 - exp(-nu * timestep))
    immunity[k] ~ binomial(T[k], 1 - exp(-k_erlang * tau * timestep))
    loss_immunity ~ binomial(immunity[k_erlang - 1], 1 - alpha)
    S <- S - infection + loss_immunity
    E <- E + infection - incubation
    I <- I + incubation - loss_infectiousness
    T[k] <- T[k] + (k == 0 ? loss_infectiousness : 0) + (k > 0 ? immunity[k - 1] : 0) - immunity[k]
    L <- L + immunity[k_erlang - 1] - loss_immunity
    Inc <- Inc + infection
  }
  sub observation {
    Cases ~ poisson(rho * Inc)
  }
}
