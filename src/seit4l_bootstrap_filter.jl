#=
Simple SEIT4L bootstrap particle filter (no SSMProblems dependency).
Used by chain generation scripts.
=#

using Distributions
using StatsBase  # For wsample

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
