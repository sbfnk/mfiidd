using Random
using GeneralisedFilters
using ForwardDiff

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
    # Strip ForwardDiff Duals before constructing the SSM components: the
    # bootstrap filter is non-differentiable, and SEIT4LDynamics /
    # PoissonObservation declare concrete Float64 fields. Mirrors the pattern
    # used in sessions/pmcmc.qmd.
    θ_f64 = Dict{Symbol,Float64}(k => ForwardDiff.value(v) for (k, v) in θ)

    # Define SSM components
    initial = SEIT4LInitial(init_state)
    dynamics = SEIT4LDynamics(θ_f64)
    observation = PoissonObservation(θ_f64[:ρ])

    # Create state-space model
    model = StateSpaceModel(initial, dynamics, observation)

    # Run bootstrap particle filter
    rng = Random.default_rng()
    algo = BF(n_particles)  # Bootstrap Filter
    _, log_lik = GeneralisedFilters.filter(rng, model, algo, obs)

    return log_lik
end
