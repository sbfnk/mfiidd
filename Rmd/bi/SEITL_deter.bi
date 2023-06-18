model seitlDeter {
  const k_erlang = 1
  const N = 1000
  dim k(k_erlang)
  state S, E, I, T[k], L, Inc
  param R_0, D_lat, D_inf, alpha, D_imm, rho
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
  sub transition {
    inline beta = R_0/D_inf
    inline epsilon = 1/D_lat
    inline nu = 1/D_inf
    inline tau = 1/D_imm
    Inc <- 0
    ode {
      dS/dt = -beta * S * I/N + (1-alpha) * tau * k_erlang * T[k_erlang - 1]
      dE/dt = beta * S * I/N - epsilon * E
      dI/dt = epsilon * E - nu * I
      dT[k]/dt = + (k == 0 ? nu * I : 0)
      - k_erlang * tau * T[k]
      + (k > 0 ? k_erlang * tau * T[k-1] : 0)
      dL/dt = alpha * k_erlang * tau * T[k_erlang - 1]
      dInc/dt = epsilon * E
    }
  }
  sub observation {
    Cases ~ poisson(rho * Inc)
  }
}
