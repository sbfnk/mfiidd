#=
SIR Model Definition
Used in sessions 1-4 (Introduction, MCMC, MCMC Diagnostics, Model Checking)

Usage:
    include(joinpath(@__DIR__, "..", "src", "sir_model.jl"))
=#

using DifferentialEquations
using Distributions
using DataFrames

"""
    sir_ode!(du, u, p, t)

SIR model ordinary differential equations.

# States
- `u[1]` (S): Susceptible
- `u[2]` (I): Infectious
- `u[3]` (R): Recovered

# Parameters
- `p[1]` (R_0): Basic reproduction number
- `p[2]` (D_inf): Infectious period (days)
"""
function sir_ode!(du, u, p, t)
    S, I, R = u
    R_0, D_inf = p
    N = S + I + R

    β = R_0 / D_inf
    ν = 1 / D_inf

    du[1] = -β * S * I / N          # dS/dt
    du[2] = β * S * I / N - ν * I   # dI/dt
    du[3] = ν * I                   # dR/dt
end

"""
    simulate_sir(R_0, D_inf, S0, I0, times)

Simulate the deterministic SIR model.

# Arguments
- `R_0`: Basic reproduction number
- `D_inf`: Infectious period (days)
- `S0`: Initial susceptible population
- `I0`: Initial infectious population
- `times`: Time points to return (e.g., 0.0:1.0:30.0)

# Returns
DataFrame with columns: time, S, I, R, Inc (daily incidence)
"""
function simulate_sir(R_0, D_inf, S0, I0, times)
    times_vec = collect(times)
    u0 = Float64[S0, I0, 0.0]
    prob = ODEProblem(sir_ode!, u0, (times_vec[1], times_vec[end]), [R_0, D_inf])
    sol = solve(prob, Tsit5(), saveat=times_vec)

    df = DataFrame(
        time = sol.t,
        S = [sol.u[i][1] for i in 1:length(sol)],
        I = [sol.u[i][2] for i in 1:length(sol)],
        R = [sol.u[i][3] for i in 1:length(sol)]
    )

    # Compute daily incidence from change in R
    df.Inc = [0.0; diff(df.R)]

    return df
end
