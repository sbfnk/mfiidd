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

# =============================================================================
# SEIT4L Particle Filter
# =============================================================================

function particle_filter_seit4l(θ, obs, n_particles; init_state=[279.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0])
    n_obs = length(obs)
    ρ = θ[:ρ]

    function rates(state)
        S, E, I, T1, T2, T3, T4, L = state
        β = θ[:R_0] / θ[:D_inf]
        ϵ = 1.0 / θ[:D_lat]
        ν = 1.0 / θ[:D_inf]
        τ = 4.0 / θ[:D_imm]
        N = S + E + I + T1 + T2 + T3 + T4 + L
        [β*S*I/N, ϵ*E, ν*I, τ*T1, τ*T2, τ*T3, (1-θ[:α])*τ*T4, θ[:α]*τ*T4]
    end

    stoich = [[-1,1,0,0,0,0,0,0], [0,-1,1,0,0,0,0,0], [0,0,-1,1,0,0,0,0],
              [0,0,0,-1,1,0,0,0], [0,0,0,0,-1,1,0,0], [0,0,0,0,0,-1,1,0],
              [1,0,0,0,0,0,-1,0], [0,0,0,0,0,0,-1,1]]

    function gillespie_step!(state, dt)
        t, incidence = 0.0, 0
        while t < dt
            r = rates(state)
            total = sum(r)
            total ≤ 0 && break
            τ_step = randexp() / total
            t + τ_step > dt && break
            t += τ_step
            cum, rnd, ev = 0.0, rand()*total, 0
            for i in 1:8
                cum += r[i]
                if rnd ≤ cum; ev = i; break; end
            end
            for j in 1:8; state[j] += stoich[ev][j]; end
            ev == 2 && (incidence += 1)
        end
        incidence
    end

    particles = [copy(init_state) for _ in 1:n_particles]
    log_lik = 0.0

    for t in 1:n_obs
        inc = [gillespie_step!(particles[i], 1.0) for i in 1:n_particles]
        log_w = [logpdf(Poisson(max(ρ*inc[i], 1e-10)), obs[t]) for i in 1:n_particles]
        max_lw = maximum(log_w)
        w = exp.(log_w .- max_lw)
        log_lik += max_lw + log(mean(w))
        w ./= sum(w)
        ess = 1.0 / sum(w.^2)
        if ess < n_particles / 2
            idx = wsample(1:n_particles, w, n_particles)
            particles = [copy(particles[i]) for i in idx]
        end
    end
    log_lik
end

# =============================================================================
# PMMH Model
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
