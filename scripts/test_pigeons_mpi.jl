# Test Pigeons.jl with MPI parallelization
# Run with: julia --project=. scripts/test_pigeons_mpi.jl

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
# Particle filter for SEIT4L (fewer particles for speed)
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
# Run with MPI (local processes for testing)
# =============================================================================

n_particles = 32  # Fewer particles for speed
n_rounds = 8      # More rounds for more samples

println("Creating model with $n_particles particles...")
model = pmmh_seit4l(flu_tdc.obs, n_particles)

println("Running Pigeons with 4 local MPI processes...")
println("This distributes the parallel tempering chains across processes.")

pt = pigeons(
    target = TuringLogPotential(model),
    n_rounds = n_rounds,
    on = ChildProcess(
        n_local_mpi_processes = 4,
        n_threads = 1
    ),
    record = [traces]
)

println("\nPigeons result:")
println(pt)

# Get samples
samples = Chains(pt)
println("\nSamples summary:")
println(samples)

# Save results
using JLD2
output_path = datadir("pigeons_test_chain.jld2")
@save output_path samples
println("\nSaved to $output_path")
