using Random
using GeneralisedFilters

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
