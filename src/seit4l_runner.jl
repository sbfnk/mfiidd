using Random: default_rng
using GeneralisedFilters: BF, filter
using ForwardDiff: value
using SSMProblems: StateSpaceModel

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
    # Coerce init_state to Float64 so callers can pass an integer vector
    # (e.g. [279, 0, 2, ...]) without hitting SEIT4LInitial's Vector{Float64}
    # signature.
    init_state_f64 = collect(Float64.(init_state))
    length(init_state_f64) == 8 ||
        throw(ArgumentError("init_state must be [S, E, I, T1, T2, T3, T4, L]"))

    # Strip ForwardDiff Duals before constructing the SSM components: the
    # bootstrap filter is non-differentiable, and SEIT4LDynamics /
    # PoissonObservation declare concrete Float64 fields. Mirrors the pattern
    # used in sessions/pmcmc.qmd.
    θ_f64 = Dict{Symbol,Float64}(k => value(v) for (k, v) in θ)

    # Define SSM components
    initial = SEIT4LInitial(init_state_f64)
    dynamics = SEIT4LDynamics(θ_f64)
    observation = PoissonObservation(θ_f64[:ρ])

    # Create state-space model
    model = StateSpaceModel(initial, dynamics, observation)

    # Run bootstrap particle filter
    rng = default_rng()
    algo = BF(n_particles)  # Bootstrap Filter
    _, log_lik = filter(rng, model, algo, obs)

    return log_lik
end
