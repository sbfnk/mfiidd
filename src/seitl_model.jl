#=
SEITL and SEIT4L Model Definitions
Used in sessions 5+ (SEITL, Particle Filters, PMCMC, etc.)

Models:
- SEITL: SIR with Exposed compartment and Temporary/Long-term immunity
- SEIT4L: SEITL with Erlang-distributed temporary immunity (4 sub-stages)

Usage:
    include(joinpath(@__DIR__, "..", "src", "seitl_model.jl"))
=#

using DifferentialEquations
using Distributions
using DataFrames
using Random

# =============================================================================
# SEITL Model
# =============================================================================

"""
    seitl_ode!(du, u, p, t)

SEITL model: SIR with Exposed compartment and Temporary/Long-term immunity.

# States
- S: Susceptible
- E: Exposed (latent)
- I: Infectious
- T: Temporary immunity
- L: Long-term immunity
- Inc: Cumulative incidence (for tracking)

# Parameters
- R_0: Basic reproduction number
- D_lat: Latent period (days)
- D_inf: Infectious period (days)
- α: Probability of developing long-term immunity
- D_imm: Duration of temporary immunity (days)
"""
function seitl_ode!(du, u, p, t)
    R_0, D_lat, D_inf, α, D_imm = p

    β = R_0 / D_inf
    ϵ = 1.0 / D_lat
    ν = 1.0 / D_inf
    τ = 1.0 / D_imm

    S, E, I, T, L, Inc = u
    N = S + E + I + T + L

    du[1] = -β * S * I / N + (1 - α) * τ * T  # dS/dt
    du[2] = β * S * I / N - ϵ * E              # dE/dt
    du[3] = ϵ * E - ν * I                      # dI/dt
    du[4] = ν * I - τ * T                      # dT/dt
    du[5] = α * τ * T                          # dL/dt
    du[6] = ϵ * E                              # dInc/dt (cumulative)
end

"""
    simulate_seitl_deterministic(θ, init_state, times)

Simulate the deterministic SEITL model.

# Arguments
- `θ`: Dict with keys :R_0, :D_lat, :D_inf, :α, :D_imm
- `init_state`: Dict with keys :S, :E, :I, :T, :L
- `times`: Time points to return

# Returns
DataFrame with columns: time, S, E, I, T, L, Inc (daily incidence)
"""
function simulate_seitl_deterministic(θ, init_state, times)
    params = [θ[:R_0], θ[:D_lat], θ[:D_inf], θ[:α], θ[:D_imm]]

    u0 = [init_state[:S], init_state[:E], init_state[:I],
          init_state[:T], init_state[:L], 0.0]

    prob = ODEProblem(seitl_ode!, u0, (times[1], times[end]), params)
    sol = solve(prob, Tsit5(), saveat=times)

    state_names = [:S, :E, :I, :T, :L, :Inc_cumulative]
    sol_matrix = reduce(hcat, sol.u)'
    df = DataFrame(sol_matrix, state_names)
    df[!, :time] = sol.t

    # Convert cumulative to daily incidence
    df[!, :Inc] = [0.0; diff(df.Inc_cumulative)]

    return df
end

"""
    simulate_seitl_stochastic(θ, init_state, times)

Simulate the stochastic SEITL model using the Gillespie algorithm.

# Arguments
- `θ`: Dict with keys :R_0, :D_lat, :D_inf, :α, :D_imm
- `init_state`: Dict with keys :S, :E, :I, :T, :L
- `times`: Time points (assumed daily: 0, 1, 2, ...)

# Returns
DataFrame with columns: time, S, E, I, T, L, Inc (daily incidence)
"""
function simulate_seitl_stochastic(θ, init_state, times)
    R_0, D_lat, D_inf, α, D_imm = θ[:R_0], θ[:D_lat], θ[:D_inf], θ[:α], θ[:D_imm]
    β = R_0 / D_inf
    ϵ = 1.0 / D_lat
    ν = 1.0 / D_inf
    τ = 1.0 / D_imm

    # State: [S, E, I, T, L]
    state = Float64[init_state[:S], init_state[:E], init_state[:I],
                    init_state[:T], init_state[:L]]

    # Stoichiometry: how each transition changes [S, E, I, T, L]
    stoich = [
        [-1, 1, 0, 0, 0],   # S → E
        [0, -1, 1, 0, 0],   # E → I
        [0, 0, -1, 1, 0],   # I → T
        [1, 0, 0, -1, 0],   # T → S
        [0, 0, 0, -1, 1]    # T → L
    ]

    function rates(s)
        S, E, I, T, L = s
        N = S + E + I + T + L
        [β*S*I/N, ϵ*E, ν*I, (1-α)*τ*T, α*τ*T]
    end

    # Gillespie step: simulate one time unit, return daily incidence
    function gillespie_step!(s, dt)
        t, daily_inc = 0.0, 0
        while t < dt
            r = rates(s)
            total_rate = sum(r)
            total_rate ≤ 0 && break
            τ_wait = randexp() / total_rate
            t + τ_wait > dt && break
            t += τ_wait
            # Select event
            cum, rnd, event = 0.0, rand() * total_rate, 0
            for i in 1:5
                cum += r[i]
                if rnd ≤ cum
                    event = i
                    break
                end
            end
            # Apply transition
            for j in 1:5
                s[j] += stoich[event][j]
            end
            # Track E→I transitions as new cases
            event == 2 && (daily_inc += 1)
        end
        daily_inc
    end

    # Simulate day-by-day
    n_days = length(times)
    results = DataFrame(time=collect(times), S=zeros(n_days), E=zeros(n_days),
                       I=zeros(n_days), T=zeros(n_days), L=zeros(n_days),
                       Inc=zeros(n_days))

    for (i, t) in enumerate(times)
        results.S[i], results.E[i], results.I[i] = state[1], state[2], state[3]
        results.T[i], results.L[i] = state[4], state[5]
        if i < n_days
            results.Inc[i+1] = gillespie_step!(state, times[i+1] - t)
        end
    end

    return results
end

# =============================================================================
# SEIT4L Model (Erlang-distributed temporary immunity)
# =============================================================================

"""
    seit4l_ode!(du, u, p, t)

SEIT4L model: SEITL with T compartment split into 4 stages (Erlang-4 distribution).

# States
S, E, I, T1, T2, T3, T4, L, Inc (cumulative)

# Parameters
Same as SEITL: R_0, D_lat, D_inf, α, D_imm
"""
function seit4l_ode!(du, u, p, t)
    R_0, D_lat, D_inf, α, D_imm = p

    β = R_0 / D_inf
    ϵ = 1.0 / D_lat
    ν = 1.0 / D_inf
    τ = 4.0 / D_imm  # Rate through each T sub-stage

    S, E, I, T1, T2, T3, T4, L, Inc = u
    N = S + E + I + T1 + T2 + T3 + T4 + L

    du[1] = -β * S * I / N + (1 - α) * τ * T4  # dS/dt
    du[2] = β * S * I / N - ϵ * E              # dE/dt
    du[3] = ϵ * E - ν * I                      # dI/dt
    du[4] = ν * I - τ * T1                     # dT1/dt
    du[5] = τ * T1 - τ * T2                    # dT2/dt
    du[6] = τ * T2 - τ * T3                    # dT3/dt
    du[7] = τ * T3 - τ * T4                    # dT4/dt
    du[8] = α * τ * T4                         # dL/dt
    du[9] = ϵ * E                              # dInc/dt (cumulative)
end

"""
    simulate_seit4l_deterministic(θ, init_state, times)

Simulate the deterministic SEIT4L model.

# Arguments
- `θ`: Dict with keys :R_0, :D_lat, :D_inf, :α, :D_imm
- `init_state`: Dict with keys :S, :E, :I, :T1, :T2, :T3, :T4, :L
- `times`: Time points to return

# Returns
DataFrame with columns: time, S, E, I, T1, T2, T3, T4, L, Inc
"""
function simulate_seit4l_deterministic(θ, init_state, times)
    params = [θ[:R_0], θ[:D_lat], θ[:D_inf], θ[:α], θ[:D_imm]]

    u0 = [init_state[:S], init_state[:E], init_state[:I],
          init_state[:T1], init_state[:T2], init_state[:T3], init_state[:T4],
          init_state[:L], 0.0]

    prob = ODEProblem(seit4l_ode!, u0, (times[1], times[end]), params)
    sol = solve(prob, Tsit5(), saveat=times)

    state_names = [:S, :E, :I, :T1, :T2, :T3, :T4, :L, :Inc_cumulative]
    sol_matrix = reduce(hcat, sol.u)'
    df = DataFrame(sol_matrix, state_names)
    df[!, :time] = sol.t

    df[!, :Inc] = [0.0; diff(df.Inc_cumulative)]

    return df
end

# =============================================================================
# Observation Process Helper
# =============================================================================

"""
    generate_observations(trajectory, ρ)

Generate observed incidence using Poisson observation process.

# Arguments
- `trajectory`: DataFrame with Inc column
- `ρ`: Reporting rate (0-1)

# Returns
Vector of observed counts
"""
function generate_observations(trajectory, ρ)
    return [rand(Poisson(max(ρ * inc, 1e-10))) for inc in trajectory.Inc]
end
