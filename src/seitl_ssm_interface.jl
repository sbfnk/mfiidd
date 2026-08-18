using Random
using Distributions
using SSMProblems

"""
SEITL latent dynamics.
State vector: [S, E, I, T, L, daily_inc]
The last element tracks daily incidence for the observation process.
"""
struct SEITLDynamics <: SSMProblems.LatentDynamics
    θ::Dict{Symbol, Float64}
end

function SSMProblems.simulate(rng::AbstractRNG, dyn::SEITLDynamics,
                              step::Integer, prev_state; kwargs...)
    state = collect(prev_state[1:5])  # Extract compartments
    daily_inc = gillespie_step_seitl!(rng, state, dyn.θ)
    return vcat(state, daily_inc)     # Append daily incidence
end

"""
Initial state distribution (deterministic).
"""
struct SEITLInitial <: SSMProblems.StatePrior
    init_state::Vector{Float64}
end

function SSMProblems.simulate(rng::AbstractRNG, prior::SEITLInitial; kwargs...)
    vcat(prior.init_state, 0.0)  # Append 0 for initial daily incidence
end
