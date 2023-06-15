## load Tristan da Cunha data
data(fluTdc1971)

## define deterministic skeleton
SEITL.skel.c <- '
    double trans[5];

    double beta = R0 / D_inf;
    double epsilon = 1 / D_lat;
    double nu = 1 / D_inf;
    double tau = 1 / D_imm;

    double N = S + E + I + T + L;

    trans[0] = beta * I / N * S;
    trans[1] = epsilon * E;
    trans[2] = nu * I;
    trans[3] = alpha * tau * T;
    trans[4] = (1 - alpha) * tau * T;

    DS = -trans[0] + trans[4];
    DE = trans[0] - trans[1];
    DI = trans[1] - trans[2];
    DT = trans[2] - trans[3] - trans[4];
    DL = trans[3];
    DInc = trans[1];
'

## define stochastic model, for use with euler.sim, see ?euler.sim
SEITL.sim.c <- '
    double rate[5];
    double dN[5];

    double beta = R0 / D_inf;
    double epsilon = 1 / D_lat;
    double nu = 1 / D_inf;
    double tau = 1 / D_imm;

    double N = S + E + I + T + L;

    rate[0] = beta * I / N;
    rate[1] = epsilon;
    rate[2] = nu;
    rate[3] = alpha * tau;
    rate[4] = (1 - alpha) * tau;

    reulermultinom(1, S, &rate[0], dt, &dN[0]);
    reulermultinom(1, E, &rate[1], dt, &dN[1]);
    reulermultinom(1, I, &rate[2], dt, &dN[2]);
    reulermultinom(2, T, &rate[3], dt, &dN[3]);

    S += -dN[0] + dN[4];
    E += dN[0] - dN[1];
    I += dN[1] - dN[2];
    T += dN[2] - dN[3] - dN[4];
    L += dN[3];
    Inc += dN[1];
'

## define sampling random point observations
SEITL.rmeas.c <- '
    obs = rpois(rho * Inc > 0 ? rho * Inc : 0);
'

## define point observation probability density
SEITL.dmeas.c <- '
    lik = dpois(obs, rho * Inc > 0 ? rho * Inc : 0, give_log);
'

## define prior density
SEITL.dprior.c <- '
  lik = dunif(R0, 1, 50, 1) +
          dunif(D_lat, 0, 10, 1) +
          dunif(D_inf, 0, 15, 1) +
          dunif(D_imm, 0, 50, 1) +
          dunif(alpha, 0, 1, 1) +
          dunif(rho, 0, 1, 1);

  lik = give_log ? lik : exp(lik);
'

SEITL.logtrans.c <- '
    TR0 = log(R0);
    TD_inf = log(D_inf);
    TD_lat = log(D_lat);
    TD_imm = log(D_imm);
    Talpha = log(alpha);
    Trho = log(rho);
'

SEITL.exptrans.c <- '
    TR0 = exp(R0);
    TD_inf = exp(D_inf);
    TD_lat = exp(D_lat);
    TD_imm = exp(D_imm);
    Talpha = exp(alpha);
    Trho = exp(rho);
'

## construct pomp object
SEITL_pomp <- pomp(data = fluTdc1971[, c("time", "obs")],
                   skeleton = vectorfield(Csnippet(SEITL.skel.c)),
                   rprocess = euler.sim(step.fun = Csnippet(SEITL.sim.c),
                                        delta.t = 0.1),
                   rmeasure = Csnippet(SEITL.rmeas.c),
                   dmeasure = Csnippet(SEITL.dmeas.c),
                   dprior = Csnippet(SEITL.dprior.c), 
                   toEstimationScale = Csnippet(SEITL.logtrans.c),
                   fromEstimationScale = Csnippet(SEITL.exptrans.c), 
                   times = "time",
                   t0 = 1,
                   zeronames = "Inc", 
                   paramnames = c("R0", "D_inf", "D_lat", "D_imm", "alpha", "rho"),
                   statenames = c("S", "E", "I", "T", "L", "Inc"),
                   obsnames = c("obs"))
