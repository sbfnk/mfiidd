#=
SEIT4L Particle Filter Components

This file contains the shared Gillespie simulator and particle filter
for the SEIT4L model.

Used by:
- Particle Filters session (sessions/particle_filters.qmd)
- Particle MCMC session (sessions/pmcmc.qmd)
- Chain generation scripts (scripts/generate_pmcmc_*.jl)
=#

using Random
using Distributions
using StatsBase  # For wsample in simple bootstrap filter

"""
    gillespie_step(rng, state, θ)

Simulate SEIT4L for one day using the Gillespie algorithm.

# Arguments
- `rng`: Random number generator
- `state`: Vector [S, E, I, T1, T2, T3, T4, L]
- `θ`: Parameter dictionary with keys :R_0, :D_lat, :D_inf, :α, :D_imm

# Returns
- `(new_state, daily_incidence)`: Updated state and number of new cases
"""
function gillespie_step(rng::AbstractRNG, state::Vector{Float64}, θ::Dict)
    β = θ[:R_0] / θ[:D_inf]
    ϵ = 1.0 / θ[:D_lat]
    ν = 1.0 / θ[:D_inf]
    τ = 4.0 / θ[:D_imm]
    α = θ[:α]

    # Copy state for modification
    s = copy(state)

    # Stoichiometry: how each transition changes [S, E, I, T1, T2, T3, T4, L]
    stoich = [
        [-1, 1, 0, 0, 0, 0, 0, 0],   # S → E (infection)
        [0, -1, 1, 0, 0, 0, 0, 0],   # E → I (becoming infectious)
        [0, 0, -1, 1, 0, 0, 0, 0],   # I → T1 (recovery)
        [0, 0, 0, -1, 1, 0, 0, 0],   # T1 → T2
        [0, 0, 0, 0, -1, 1, 0, 0],   # T2 → T3
        [0, 0, 0, 0, 0, -1, 1, 0],   # T3 → T4
        [1, 0, 0, 0, 0, 0, -1, 0],   # T4 → S (immunity wanes)
        [0, 0, 0, 0, 0, 0, -1, 1]    # T4 → L (long-term immunity)
    ]

    function rates(s)
        S, E, I, T1, T2, T3, T4, L = s
        N = S + E + I + T1 + T2 + T3 + T4 + L
        [β*S*I/N, ϵ*E, ν*I, τ*T1, τ*T2, τ*T3, (1-α)*τ*T4, α*τ*T4]
    end

    # Simulate one day
    t, daily_inc = 0.0, 0
    while t < 1.0
        r = rates(s)
        total_rate = sum(r)
        total_rate ≤ 0 && break

        # Time to next event
        τ_wait = randexp(rng) / total_rate
        t + τ_wait > 1.0 && break
        t += τ_wait

        # Select which event occurs
        cum, rnd, event = 0.0, rand(rng) * total_rate, 0
        for i in 1:8
            cum += r[i]
            if rnd ≤ cum
                event = i
                break
            end
        end

        # Apply the transition
        for j in 1:8
            s[j] += stoich[event][j]
        end

        # E → I transitions count as new cases
        event == 2 && (daily_inc += 1)
    end

    return s, daily_inc
end

"""
    gillespie_step_seit4l!(state, θ, dt)

In-place version for simple bootstrap filter (no RNG argument, uses global RNG).
"""
function gillespie_step_seit4l!(state::Vector{Float64}, θ::Dict, dt::Float64=1.0)
    new_state, inc = gillespie_step(Random.default_rng(), state, θ)
    for i in 1:8
        state[i] = new_state[i]
    end
    return inc
end

"""
    particle_filter_seit4l(θ, obs, n_particles; init_state)

Simple bootstrap particle filter for SEIT4L (no SSMProblems dependency).
Used by chain generation scripts.

# Arguments
- `θ`: Parameter dictionary (must include :R_0, :D_lat, :D_inf, :α, :D_imm, :ρ)
- `obs`: Vector of observed daily incidence
- `n_particles`: Number of particles
- `init_state`: Initial state vector [S, E, I, T1, T2, T3, T4, L]

# Returns
- `log_likelihood`: Estimated log-likelihood
"""
function particle_filter_seit4l(θ, obs, n_particles;
                                init_state=[279.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0])
    n_obs = length(obs)
    ρ = θ[:ρ]

    particles = [copy(init_state) for _ in 1:n_particles]
    log_lik = 0.0

    for t in 1:n_obs
        # Propagate particles
        inc = [gillespie_step_seit4l!(particles[i], θ) for i in 1:n_particles]

        # Weight by observation likelihood
        log_w = [logpdf(Poisson(max(ρ * inc[i], 1e-10)), obs[t]) for i in 1:n_particles]

        # Log-sum-exp trick for numerical stability
        max_lw = maximum(log_w)
        w = exp.(log_w .- max_lw)
        log_lik += max_lw + log(mean(w))

        # Normalize weights
        w ./= sum(w)

        # Resample if ESS too low
        ess = 1.0 / sum(w.^2)
        if ess < n_particles / 2
            idx = wsample(1:n_particles, Weights(w), n_particles)
            particles = [copy(particles[i]) for i in idx]
        end
    end

    return log_lik
end

#=
SSMProblems.jl Interface

These types define the SEIT4L model as a state-space model for use with
GeneralisedFilters.jl particle filtering algorithms.

Note: Sessions using this interface must also load SSMProblems and GeneralisedFilters.
=#

using SSMProblems
using GeneralisedFilters

"""
SEIT4L latent dynamics.
State vector: [S, E, I, T1, T2, T3, T4, L, daily_inc]
The last element tracks daily incidence for the observation process.
"""
struct SEIT4LDynamics <: SSMProblems.LatentDynamics
    θ::Dict{Symbol, Float64}
end

function SSMProblems.simulate(rng::AbstractRNG, dyn::SEIT4LDynamics,
                              step::Integer, prev_state; kwargs...)
    # Extract compartments (first 8 elements)
    state = collect(prev_state[1:8])

    # Simulate one day
    new_state, daily_inc = gillespie_step(rng, state, dyn.θ)

    # Return state with daily incidence appended
    return vcat(new_state, daily_inc)
end

"""
Poisson observation process.
Observes daily incidence (last element of state) with reporting rate ρ.
"""
struct PoissonObservation <: SSMProblems.ObservationProcess
    ρ::Float64
end

function SSMProblems.distribution(obs::PoissonObservation, step::Integer,
                                  state; kwargs...)
    daily_inc = state[end]
    Poisson(max(obs.ρ * daily_inc, 1e-10))
end

"""
Initial state distribution (deterministic).
"""
struct SEIT4LInitial <: SSMProblems.StatePrior
    init_state::Vector{Float64}
end

function SSMProblems.simulate(rng::AbstractRNG, prior::SEIT4LInitial; kwargs...)
    vcat(prior.init_state, 0.0)  # Append 0 for initial daily incidence
end

"""
    run_particle_filter(θ, obs, n_particles; init_state)

Run bootstrap particle filter for SEIT4L and return log-likelihood.

# Arguments
- `θ`: Parameter dictionary (must include :R_0, :D_lat, :D_inf, :α, :D_imm, :ρ)
- `obs`: Vector of observed daily incidence
- `n_particles`: Number of particles
- `init_state`: Initial state vector [S, E, I, T1, T2, T3, T4, L]

# Returns
- `log_likelihood`: Estimated log-likelihood
"""
function run_particle_filter(θ, obs, n_particles;
                             init_state=[279.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0])
    # Define SSM components
    initial = SEIT4LInitial(init_state)
    dynamics = SEIT4LDynamics(θ)
    observation = PoissonObservation(θ[:ρ])

    # Create state-space model
    model = StateSpaceModel(initial, dynamics, observation)

    # Run bootstrap particle filter
    rng = Random.default_rng()
    algo = BF(n_particles)  # Bootstrap Filter
    _, log_lik = GeneralisedFilters.filter(rng, model, algo, obs)

    return log_lik
end
