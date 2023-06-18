## define deterministic skeleton
seit2lDeterSkel <- '
    double trans[6];

    double beta = R0 / D_inf;
    double epsilon = 1 / D_lat;
    double nu = 1 / D_inf;
    double tau = 1 / D_imm;

    double N = S + E + I + T1 + T2 + L;

    trans[0] = beta * I / N * S;
    trans[1] = epsilon * E;
    trans[2] = nu * I;
    trans[3] = 2 * tau * T1;
    trans[4] = 2 * alpha * tau * T2;
    trans[5] = 2 * (1 - alpha) * tau * T2;

    DS = -trans[0] + trans[5];
    DE = trans[0] - trans[1];
    DI = trans[1] - trans[2];
    DT1 = trans[2] - trans[3];
    DT2 = trans[3] - trans[4] - trans[5];
    DL = trans[4];
    DInc = trans[1];
'

## define stochastic model, for use with euler, see ?euler
seit2lStochSim <- '
    double rate[6];
    double dN[6];

    double beta = R0 / D_inf;
    double epsilon = 1 / D_lat;
    double nu = 1 / D_inf;
    double tau = 1 / D_imm;

    double N = S + E + I + T1 + T2 + L;

    rate[0] = beta * I / N;
    rate[1] = epsilon;
    rate[2] = nu;
    rate[3] = 2 * tau;
    rate[4] = 2 * alpha * tau;
    rate[5] = 2 * (1 - alpha) * tau;

    reulermultinom(1, S, &rate[0], dt, &dN[0]);
    reulermultinom(1, E, &rate[1], dt, &dN[1]);
    reulermultinom(1, I, &rate[2], dt, &dN[2]);
    reulermultinom(1, T1, &rate[3], dt, &dN[3]);
    reulermultinom(2, T2, &rate[4], dt, &dN[4]);

    S += -dN[0] + dN[5];
    E += dN[0] - dN[1];
    I += dN[1] - dN[2];
    T1 += dN[2] - dN[3];
    T2 += dN[3] - dN[4] - dN[5];
    L += dN[4];
    Inc += dN[1];
'


## construct pomp object
seit2lPomp <- pomp(
  data = fluTdc1971[, c("time", "obs")],
  skeleton = vectorfield(seit2lDeterSkel),
  rprocess = euler(step.fun = seit2lStochSim, delta.t = 0.1),
  rmeasure = seitlGenObsPoint,
  dmeasure = seitlPointLike,
  dprior = seitlPrior,
  partrans = parameter_trans(
    log = c("R0", "D_inf", "D_lat", "D_imm", "alpha", "rho")
  ),
  times = "time",
  t0 = 1,
  accumvars = "Inc",
  paramnames = c("R0", "D_inf", "D_lat", "D_imm", "alpha", "rho"),
  statenames = c("S", "E", "I", "T1", "T2", "L", "Inc"),
  obsnames = c("obs")
)
