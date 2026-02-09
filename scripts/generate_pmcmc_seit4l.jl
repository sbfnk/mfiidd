# Script to generate SEIT4L PMMH chain only
# Run with: julia --project=. scripts/generate_pmcmc_seit4l.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
using Distributions
using DataFrames
using Turing
using MCMCChains
using CSV
using DrWatson
using StatsBase
using AdvancedMH

Random.seed!(1234)

# Load data
flu_tdc = CSV.read(datadir("flu_tdc_1971.csv"), DataFrame)

# Load shared particle filter
include(joinpath(@__DIR__, "..", "src", "seit4l_particle_filter.jl"))

# =============================================================================
# PMMH Model
# =============================================================================

@model function pmmh_seit4l(obs, n_particles)
    R_0 ~ truncated(Normal(3.0, 2.0), lower=1.0)
    D_lat ~ truncated(Normal(2.0, 1.0), lower=0.5)
    D_inf ~ truncated(Normal(3.0, 2.0), lower=0.5)
    α ~ Beta(2, 2)
    D_imm ~ truncated(Normal(15.0, 10.0), lower=1.0)
    ρ ~ Beta(2, 2)

    θ = Dict(:R_0 => R_0, :D_lat => D_lat, :D_inf => D_inf,
             :α => α, :D_imm => D_imm, :ρ => ρ)

    log_lik = particle_filter_seit4l(θ, obs, n_particles)
    Turing.@addlogprob! log_lik
end

# =============================================================================
# Run PMMH with RAM
# =============================================================================

n_particles = 256
n_samples = 500_000
burnin = 50_000
thinning = 50  # → 9000 final samples

# Helper function to save chain as CSV
function save_chain_csv(chain, path)
    params = chain.name_map.parameters
    internals = chain.name_map.internals
    all_names = vcat(params, internals)
    data_matrix = chain.value.data[:, :, 1]
    df = DataFrame(data_matrix, all_names)
    CSV.write(path, df)
end

println("=" ^ 60)
println("Running PMMH for SEIT4L with RAM")
println("  Particles: $n_particles")
println("  Total samples: $n_samples")
println("  Burnin: $burnin")
println("  Thinning: $thinning")
println("  Final samples: $((n_samples - burnin) ÷ thinning)")
println("=" ^ 60)

t_start = time()
model_seit4l = pmmh_seit4l(flu_tdc.obs, n_particles)

chain_seit4l_full = sample(
    model_seit4l,
    externalsampler(AdvancedMH.RobustAdaptiveMetropolis()),
    n_samples;
    check_model=false,
    progress=true
)

t_elapsed = time() - t_start
println("\nSEIT4L sampling took $(round(t_elapsed/60, digits=1)) minutes")

# Remove burnin and thin
chain_seit4l = chain_seit4l_full[burnin+1:thinning:end]
println("After burnin and thinning: $(length(chain_seit4l)) samples")

# Save
output_path = datadir("pmcmc_seit4l_chain.csv")
println("Saving SEIT4L chain to $output_path")
save_chain_csv(chain_seit4l, output_path)

# Print summary
println("\nSEIT4L chain summary:")
println(chain_seit4l)

println("\nEffective Sample Size:")
ess_result = ess(chain_seit4l)
for p in chain_seit4l.name_map.parameters
    println("  $p: $(round(ess_result[p, :ess], digits=1))")
end

println("\n" * "=" ^ 60)
println("Done!")
println("=" ^ 60)
