# Run Pigeons.jl parallel tempering for PMMH model
# Run with: julia --project=. scripts/run_pigeons.jl
#
# Pigeons uses parallel tempering to explore multi-modal posteriors.
# Each round doubles the number of samples (2^n_rounds total).
# n_rounds=12 gives 4096 samples.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
using Distributions
using DataFrames
using Turing
using CSV
using DrWatson
using StatsBase
using Pigeons
using MCMCChains

Random.seed!(1234)

# Load data
flu_tdc = CSV.read(datadir("flu_tdc_1971.csv"), DataFrame)

# =============================================================================
# Particle filter for SEIT4L
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
# Run Pigeons
# =============================================================================

n_particles = 128  # Particles for likelihood estimation
n_rounds = 12      # 2^12 = 4096 samples

println("Creating model with $n_particles particles...")
model = pmmh_seit4l(flu_tdc.obs, n_particles)

println("Running Pigeons with $n_rounds rounds ($(2^n_rounds) samples)...")
println("This may take several hours...")

t_start = time()
pt = pigeons(
    target = TuringLogPotential(model),
    n_rounds = n_rounds,
    record = [traces]
)
t_elapsed = time() - t_start

println("\n=== Results ===")
println("Total time: $(round(t_elapsed/60, digits=1)) minutes")
println(pt)

# Get samples as MCMCChains
samples = Chains(pt)
println("\nChain summary:")
println(samples)

# ESS
println("\nEffective Sample Size:")
ess_result = ess(samples)
for p in [:R_0, :D_lat, :D_inf, :α, :D_imm, :ρ]
    println("  $p: $(round(ess_result[p, :ess], digits=1))")
end

# Save as CSV
params = [:R_0, :D_lat, :D_inf, :α, :D_imm, :ρ]
data_matrix = Array(samples[:, params, :])[:, :, 1]
df = DataFrame(data_matrix, params)
output_path = datadir("pigeons_seit4l_samples.csv")
CSV.write(output_path, df)
println("\nSaved $(nrow(df)) samples to $output_path")
