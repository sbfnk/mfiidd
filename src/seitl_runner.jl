using Random: default_rng
using GeneralisedFilters: BF
import GeneralisedFilters
using ForwardDiff: value
using SSMProblems: StateSpaceModel

"""
    run_particle_filter_seitl(θ, obs, n_particles; init_state)

Run bootstrap particle filter for SEITL and return log-likelihood.

The SEIT4L equivalent is [`run_particle_filter`](@ref); both go through
GeneralisedFilters rather than a hand-written filtering loop, so the two models
are estimated the same way and the resampling logic is maintained upstream.

# Arguments
- `θ`: Parameter dictionary (must include :R_0, :D_lat, :D_inf, :α, :D_imm, :ρ)
- `obs`: Vector of observed daily incidence
- `n_particles`: Number of particles
- `init_state`: Initial state vector [S, E, I, T, L]

# Returns
- `log_likelihood`: Estimated log-likelihood
"""
function run_particle_filter_seitl(
    θ,
    obs,
    n_particles;
    init_state = [279.0, 0.0, 2.0, 3.0, 0.0],
)
    init_state_f64 = collect(Float64.(init_state))
    length(init_state_f64) == 5 ||
        throw(ArgumentError("init_state must be [S, E, I, T, L]"))

    # Strip ForwardDiff Duals before building the SSM components: the bootstrap
    # filter is not differentiable, and the component types declare concrete
    # Float64 fields.
    θ_f64 = Dict{Symbol, Float64}(k => value(v) for (k, v) in θ)

    model = StateSpaceModel(
        SEITLInitial(init_state_f64),
        SEITLDynamics(θ_f64),
        PoissonObservation(θ_f64[:ρ]),
    )

    _, log_lik = GeneralisedFilters.filter(default_rng(), model, BF(n_particles), obs)

    return log_lik
end
