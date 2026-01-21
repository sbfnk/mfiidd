# Script to generate pre-computed PMMH chains for the PMCMC session
# Run with: julia --project=. scripts/generate_pmcmc_chains.jl

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

# =============================================================================
# Load shared particle filter implementations
# =============================================================================

include(joinpath(@__DIR__, "..", "src", "seit4l_particle_filter.jl"))
include(joinpath(@__DIR__, "..", "src", "seitl_particle_filter.jl"))

# =============================================================================
# PMMH Models
# =============================================================================

@model function pmmh_seit4l(obs, n_particles)
    R_0 ~ Uniform(1.0, 50.0)
    D_lat ~ Uniform(0.1, 10.0)
    D_inf ~ Uniform(0.1, 15.0)
    α ~ Uniform(0.0, 1.0)
    D_imm ~ Uniform(1.0, 50.0)
    ρ ~ Uniform(0.0, 1.0)

    θ = Dict(:R_0 => R_0, :D_lat => D_lat, :D_inf => D_inf,
             :α => α, :D_imm => D_imm, :ρ => ρ)

    log_lik = particle_filter_seit4l(θ, obs, n_particles)
    Turing.@addlogprob! log_lik
end

@model function pmmh_seitl(obs, n_particles)
    R_0 ~ Uniform(1.0, 50.0)
    D_lat ~ Uniform(0.1, 10.0)
    D_inf ~ Uniform(0.1, 15.0)
    α ~ Uniform(0.0, 1.0)
    D_imm ~ Uniform(1.0, 50.0)
    ρ ~ Uniform(0.0, 1.0)

    θ = Dict(:R_0 => R_0, :D_lat => D_lat, :D_inf => D_inf,
             :α => α, :D_imm => D_imm, :ρ => ρ)

    log_lik = particle_filter_seitl(θ, obs, n_particles)
    Turing.@addlogprob! log_lik
end

# =============================================================================
# Run PMMH with Robust Adaptive Metropolis (RAM)
# =============================================================================

n_particles = 256
n_samples = 500_000  # Very long chain - RAM needs time to adapt
burnin = 50_000      # Discard first 50k as burnin
thinning = 50        # Keep every 50th sample → 9000 final samples

# Helper function to save chain as CSV
function save_chain_csv(chain, path)
    params = chain.name_map.parameters
    internals = chain.name_map.internals
    all_names = vcat(params, internals)
    data_matrix = chain.value.data[:, :, 1]
    df = DataFrame(data_matrix, all_names)
    CSV.write(path, df)
end

# Helper function to print diagnostics
function print_diagnostics(chain, name)
    println("\n$name chain summary:")
    println(chain)

    # ESS
    println("\nEffective Sample Size:")
    ess_result = ess(chain)
    for p in chain.name_map.parameters
        println("  $p: $(round(ess_result[p, :ess], digits=1))")
    end

    # Acceptance rate - use the logjoint column
    if :logjoint in chain.name_map.internals
        lj = vec(chain[:logjoint].data)
        n_accepts = sum(lj[2:end] .!= lj[1:end-1])
        accept_rate = n_accepts / (length(lj) - 1)
        println("\nAcceptance rate: $(round(accept_rate * 100, digits=1))%")
    end
end

# =============================================================================
# SEIT4L
# =============================================================================
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

# Use RAM - automatically adapts proposal covariance
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

print_diagnostics(chain_seit4l, "SEIT4L")

# =============================================================================
# SEITL
# =============================================================================
println("\n" * "=" ^ 60)
println("Running PMMH for SEITL with RAM")
println("  Particles: $n_particles")
println("  Total samples: $n_samples")
println("  Burnin: $burnin")
println("  Thinning: $thinning")
println("  Final samples: $((n_samples - burnin) ÷ thinning)")
println("=" ^ 60)

t_start = time()
model_seitl = pmmh_seitl(flu_tdc.obs, n_particles)

chain_seitl_full = sample(
    model_seitl,
    externalsampler(AdvancedMH.RobustAdaptiveMetropolis()),
    n_samples;
    check_model=false,
    progress=true
)

t_elapsed = time() - t_start
println("\nSEITL sampling took $(round(t_elapsed/60, digits=1)) minutes")

# Remove burnin and thin
chain_seitl = chain_seitl_full[burnin+1:thinning:end]
println("After burnin and thinning: $(length(chain_seitl)) samples")

# Save
output_path_seitl = datadir("pmcmc_seitl_chain.csv")
println("Saving SEITL chain to $output_path_seitl")
save_chain_csv(chain_seitl, output_path_seitl)

print_diagnostics(chain_seitl, "SEITL")

println("\n" * "=" ^ 60)
println("All done!")
println("=" ^ 60)
