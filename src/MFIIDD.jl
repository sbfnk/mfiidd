module MFIIDD

using DataFrames
using DifferentialEquations
using Distributions
using GeneralisedFilters
using PrecompileTools
using Random
using SSMProblems
using StatsBase

include("sir_model.jl")
include("seitl_model.jl")
include("seitl_particle_filter.jl")
include("seit4l_particle_filter.jl")

export sir_ode!, simulate_sir
export seitl_ode!, simulate_seitl_deterministic, simulate_seitl_stochastic
export seit4l_ode!, simulate_seit4l_deterministic
export generate_observations
export gillespie_step, gillespie_step_seitl!, gillespie_step_seit4l!
export particle_filter_seitl, particle_filter_seit4l, run_particle_filter
export SEIT4LDynamics, SEIT4LInitial, PoissonObservation

@compile_workload begin
    sir_df = simulate_sir(2.0, 4.0, 999.0, 1.0, 0.0:1.0:10.0)

    θ_seitl = Dict(:R_0 => 2.0, :D_lat => 2.0, :D_inf => 4.0,
                   :α => 0.5, :D_imm => 10.0, :ρ => 0.5)
    init_seitl = Dict(:S => 99.0, :E => 0.0, :I => 1.0, :T => 0.0, :L => 0.0)
    seitl_df = simulate_seitl_deterministic(θ_seitl, init_seitl, 0.0:1.0:5.0)
    obs = generate_observations(seitl_df, 0.5)

    init_seit4l = [99.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    particle_filter_seit4l(θ_seitl, obs, 20; init_state=init_seit4l)
    run_particle_filter(θ_seitl, obs, 20; init_state=init_seit4l)
end

end # module MFIIDD
