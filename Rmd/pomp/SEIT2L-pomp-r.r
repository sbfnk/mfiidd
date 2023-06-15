SEIT2L.skel.r <- function(x, t, params, ...) {

    with(as.list(c(x, params)),
    {
        beta <- R0 / D_inf;
        epsilon <- 1 / D_lat;
        nu <- 1 / D_inf;
        tau <- 1 / D_imm;

        N <- S + E + I + T1 + T2 + L;

        trans <- c(beta * I / N * S,
                   epsilon * E,
                   nu * I,
                   2 * tau * T1,
                   2 * alpha * tau * T2,
                   2 * (1 - alpha) * tau * T2)

        DS <- -trans[1] + trans[6];
        DE <- trans[1] - trans[2];
        DI <- trans[2] - trans[3];
        DT1 <- trans[3] - trans[4];
        DT2 <- trans[4] - trans[5] - trans[6];
        DL <- trans[5];
        DInc <- trans[2];

        c(DS, DE, DI, DT1, DT2, DL, DInc)
    })
}

SEIT2L.sim.r <- function(x, t, params, delta.t, ...)
{
    with(as.list(c(x, params)),
    {
        beta = R0 / D_inf;
        epsilon = 1 / D_lat;
        nu = 1 / D_inf;
        tau = 1 / D_imm;

        N = S + E + I + T1 + T2 + L;

        rate <- c(beta * I / N,
                  epsilon,
                  nu,
                  2 * tau * T1,
                  2 * alpha * tau * T2,
                  2 * (1 - alpha) * tau * T2)

        dN <- c(reulermultinom(n = 1, size = S, rate = rate[1], dt = delta.t),
                reulermultinom(n = 1, size = E, rate = rate[2], dt = delta.t),
                reulermultinom(n = 1, size = I, rate = rate[3], dt = delta.t),
                reulermultinom(n = 1, size = T1, rate = rate[4], dt = delta.t),
                reulermultinom(n = 2, size = T2, rate = rate[5:6], dt = delta.t))

        S <- S - dN[1] + dN[6]
        E <- E + dN[1] - dN[2]
        I <- I + dN[2] - dN[3]
        T1 <- T1 + dN[3] - dN[4]
        T2 <- T2 + dN[4] - dN[5] - dN[6];
        L <- L + dN[5];
        Inc <- Inc + dN[2]

        c(S, E, I, T1, T2, L, Inc)
     })
}

## construct pomp object
SEIT2L_pomp_r <- pomp(data = fluTdc1971[, c("time", "obs")],
                      skeleton = vectorfield(SEIT2L.skel.r),
                      rprocess = euler.sim(step.fun = SEIT2L.sim.r,
                                           delta.t = 0.1),
                      rmeasure = SEITL.rmeas.r,
                      dmeasure = SEITL.dmeas.r,
                      dprior = SEITL.dprior.r, 
                      toEstimationScale = SEITL.logtrans.r,
                      fromEstimationScale = SEITL.exptrans.r, 
                      times = "time",
                      t0 = 1,
                      zeronames = "Inc", 
                      paramnames = c("R0", "D_inf", "D_lat", "D_imm", "alpha", "rho"),
                      statenames = c("S", "E", "I", "T1", "T2", "L", "Inc"),
                      obsnames = c("obs"))

