using Random: default_rng
using GeneralisedFilters: BF, filter, DenseAncestorCallback, get_ancestry
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
function run_particle_filter(
    θ,
    obs,
    n_particles;
    init_state = [279.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0],
)
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
    θ_f64 = Dict{Symbol, Float64}(k => value(v) for (k, v) in θ)

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

"""
    filtered_incidence(θ, obs, n_particles; init_state)

Run the SEIT4L bootstrap filter and return one draw from the smoothing
distribution of daily incidence, that is a latent path conditioned on `obs`.

`DenseAncestorCallback` records the particles and their ancestor indices at
every step, so a particle drawn from the final weights can be traced back to
give its whole path. The last element of the SEIT4L state is the daily
incidence, so reading that element off the path gives the trajectory directly.

# Returns
- `Vector{Float64}` of length `length(obs)`, the filtered daily incidence
"""
function filtered_incidence(
    θ,
    obs,
    n_particles;
    init_state = [279.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0],
)
    θ_f64 = Dict{Symbol, Float64}(k => value(v) for (k, v) in θ)
    model = StateSpaceModel(
        SEIT4LInitial(collect(Float64.(init_state))),
        SEIT4LDynamics(θ_f64),
        PoissonObservation(θ_f64[:ρ]),
    )

    callback = DenseAncestorCallback(nothing)
    final, _ = filter(default_rng(), model, BF(n_particles), obs; callback)

    ## draw one particle in proportion to its final weight, then follow its
    ## ancestry back to the start
    log_w = getfield.(final.particles, :log_w)
    w = exp.(log_w .- maximum(log_w))
    u = rand() * sum(w)
    idx = findfirst(>=(u), cumsum(w))
    path = get_ancestry(callback.container, idx)

    return [path[t][end] for t in 1:length(obs)]
end
