#=
SEITL Particle Filter Components

This file contains the Gillespie simulator and bootstrap particle filter
for the SEITL model. Used by chain generation scripts.

For SEIT4L, see seit4l_particle_filter.jl
=#

using Random
using Distributions
using StatsBase

"""
    gillespie_step_seitl!(state, θ, dt)

Simulate SEITL for one time unit using the Gillespie algorithm.

# Arguments
- `state`: Vector [S, E, I, T, L] (modified in place)
- `θ`: Parameter dictionary with keys :R_0, :D_lat, :D_inf, :α, :D_imm
- `dt`: Time step (typically 1.0 for daily)

# Returns
- `incidence`: Number of new cases (E→I transitions)
"""
function gillespie_step_seitl!(state::Vector{Float64}, θ::Dict, dt::Float64=1.0)
    β = θ[:R_0] / θ[:D_inf]
    ϵ = 1.0 / θ[:D_lat]
    ν = 1.0 / θ[:D_inf]
    τ = 1.0 / θ[:D_imm]
    α = θ[:α]

    # Stoichiometry: [S, E, I, T, L]
    stoich = [
        [-1, 1, 0, 0, 0],   # S → E (infection)
        [0, -1, 1, 0, 0],   # E → I (becoming infectious)
        [0, 0, -1, 1, 0],   # I → T (recovery)
        [1, 0, 0, -1, 0],   # T → S (immunity wanes)
        [0, 0, 0, -1, 1]    # T → L (long-term immunity)
    ]

    function rates(s)
        S, E, I, T, L = s
        N = S + E + I + T + L
        [β*S*I/N, ϵ*E, ν*I, (1-α)*τ*T, α*τ*T]
    end

    t, incidence = 0.0, 0
    while t < dt
        r = rates(state)
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
            state[j] += stoich[event][j]
        end

        # E→I transitions count as new cases
        event == 2 && (incidence += 1)
    end

    return incidence
end

"""
    particle_filter_seitl(θ, obs, n_particles; init_state)

Run bootstrap particle filter for SEITL and return log-likelihood.

# Arguments
- `θ`: Parameter dictionary (must include :R_0, :D_lat, :D_inf, :α, :D_imm, :ρ)
- `obs`: Vector of observed daily incidence
- `n_particles`: Number of particles
- `init_state`: Initial state vector [S, E, I, T, L]

# Returns
- `log_likelihood`: Estimated log-likelihood
"""
function particle_filter_seitl(θ, obs, n_particles;
                               init_state=[279.0, 0.0, 2.0, 3.0, 0.0])
    n_obs = length(obs)
    ρ = θ[:ρ]

    particles = [copy(init_state) for _ in 1:n_particles]
    log_lik = 0.0

    for t in 1:n_obs
        # Propagate particles
        inc = [gillespie_step_seitl!(particles[i], θ) for i in 1:n_particles]

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
