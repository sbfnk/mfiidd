module MFIIDD

using CodeTracking: definition
using DataFrames: DataFrame
using DifferentialEquations: ODEProblem, Tsit5, solve
using Distributions: Distribution, Poisson, logpdf
using GeneralisedFilters: BF
using PrecompileTools: @compile_workload
using Random: AbstractRNG, randexp
using SSMProblems: SSMProblems, StateSpaceModel
using StatsBase: Weights, wsample

include("sir_model.jl")
include("seitl_model.jl")
include("seitl_particle_filter.jl")
include("seitl_ssm_interface.jl")
include("seitl_runner.jl")
include("seit4l_gillespie.jl")
include("seit4l_ssm_interface.jl")
include("seit4l_runner.jl")

export sir_ode!, simulate_sir
export seitl_ode!, simulate_seitl_deterministic, simulate_seitl_stochastic
export seit4l_ode!, simulate_seit4l_deterministic
export generate_observations
export gillespie_step, gillespie_step_seitl!, gillespie_step_seit4l!
export run_particle_filter, run_particle_filter_seitl, filtered_incidence
export SEIT4LDynamics, SEIT4LInitial, PoissonObservation
export SEITLDynamics, SEITLInitial
export source_for

"""
    source_for(file::AbstractString) -> String

Return the source of `src/<file>` inside the installed `MFIIDD` package.
Used by the session notebooks in place of static `{.julia include="../src/X.jl"}`
blocks so the displayed code tracks the live package rather than a path
relative to the qmd.

    source_for(f::Function) -> String

Return the source of every method of `f`, looked up via CodeTracking. Useful
when students want to inspect a single function after `using MFIIDD`.
"""
function source_for(file::AbstractString)
    base = pkgdir(@__MODULE__)
    base === nothing && error("MFIIDD package directory not found")
    path = joinpath(base, "src", file)
    isfile(path) || error("source file not found: $path")
    return read(path, String)
end

function source_for(f::Function)
    chunks = String[]
    for m in methods(f)
        d = definition(String, m)
        d === nothing && continue
        src = first(d)
        ## A default argument gives a function several methods that share one
        ## definition, so show each distinct definition once
        src in chunks || push!(chunks, src)
    end
    return join(chunks, "\n\n")
end

@compile_workload begin
    sir_df = simulate_sir(2.0, 4.0, 999.0, 1.0, 0.0:1.0:10.0)

    θ_seitl = Dict(
        :R_0 => 2.0,
        :D_lat => 2.0,
        :D_inf => 4.0,
        :α => 0.5,
        :D_imm => 10.0,
        :ρ => 0.5,
    )
    init_seitl = Dict(:S => 99.0, :E => 0.0, :I => 1.0, :T => 0.0, :L => 0.0)
    seitl_df = simulate_seitl_deterministic(θ_seitl, init_seitl, 0.0:1.0:5.0)
    obs = generate_observations(seitl_df, 0.5)

    init_seit4l = [99.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    run_particle_filter(θ_seitl, obs, 20; init_state = init_seit4l)
    run_particle_filter_seitl(θ_seitl, obs, 20; init_state = [99.0, 0.0, 1.0, 0.0, 0.0])
end

end # module MFIIDD
