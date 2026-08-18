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
function particle_filter_seit4l(θ, obs, n_particles::Integer;
                                init_state=[279.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0])
    n_particles > 0 || throw(ArgumentError("n_particles must be > 0"))
    n_obs = length(obs)
    ρ = θ[:ρ]

    particles = [copy(init_state) for _ in 1:n_particles]
    log_lik = 0.0

    # Normalised log weights, carried between steps. Resampling is only done when
    # the effective sample size drops, so on the steps in between the particles
    # are *not* equally weighted and the weights have to be remembered. Resetting
    # them each step, which is easy to do by accident, silently biases the
    # likelihood estimate and so biases whatever PMMH does with it.
    log_weights = fill(-log(n_particles), n_particles)

    for t in 1:n_obs
        # Propagate particles
        inc = [gillespie_step_seit4l!(particles[i], θ) for i in 1:n_particles]

        # Weight by observation likelihood
        log_w = [logpdf(Poisson(max(ρ * inc[i], 1e-10)), obs[t]) for i in 1:n_particles]

        # Incremental likelihood is the weighted average of the new weights, using
        # the log-sum-exp trick for numerical stability
        combined = log_weights .+ log_w
        max_lw = maximum(combined)
        increment = max_lw + log(sum(exp.(combined .- max_lw)))
        log_lik += increment

        # Renormalise
        log_weights = combined .- increment
        w = exp.(log_weights)

        # Resample if ESS too low, which resets the weights to uniform
        ess = 1.0 / sum(abs2, w)
        if ess < n_particles / 2
            idx = wsample(1:n_particles, Weights(w), n_particles)
            particles = [copy(particles[i]) for i in idx]
            log_weights = fill(-log(n_particles), n_particles)
        end
    end

    return log_lik
end
