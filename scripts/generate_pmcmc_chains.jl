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
using JLD2
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
# SEITL Particle Filter
# =============================================================================

function particle_filter_seitl(θ, obs, n_particles; init_state=[279.0, 0.0, 2.0, 3.0, 0.0])
    n_obs = length(obs)
    ρ = θ[:ρ]

    function rates(state)
        S, E, I, T, L = state
        β = θ[:R_0] / θ[:D_inf]
        ϵ = 1.0 / θ[:D_lat]
        ν = 1.0 / θ[:D_inf]
        τ = 1.0 / θ[:D_imm]
        N = S + E + I + T + L
        [β*S*I/N, ϵ*E, ν*I, (1-θ[:α])*τ*T, θ[:α]*τ*T]
    end

    stoich = [[-1,1,0,0,0], [0,-1,1,0,0], [0,0,-1,1,0], [1,0,0,-1,0], [0,0,0,-1,1]]

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
            for i in 1:5
                cum += r[i]
                if rnd ≤ cum; ev = i; break; end
            end
            for j in 1:5; state[j] += stoich[ev][j]; end
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
# Run PMMH
# =============================================================================

n_particles = 256
n_samples = 50000  # Long chain for good mixing
thinning = 10      # Keep every 10th sample

# Random walk proposals - smaller steps for ~15-25% acceptance
proposal = (
    :R_0 => AdvancedMH.RandomWalkProposal(Normal(0, 0.4)),
    :D_lat => AdvancedMH.RandomWalkProposal(Normal(0, 0.05)),
    :D_inf => AdvancedMH.RandomWalkProposal(Normal(0, 0.1)),
    :α => AdvancedMH.RandomWalkProposal(Normal(0, 0.01)),
    :D_imm => AdvancedMH.RandomWalkProposal(Normal(0, 0.4)),
    :ρ => AdvancedMH.RandomWalkProposal(Normal(0, 0.01))
)

# Initial values near posterior mode (from deterministic fit) to reduce burn-in
init_params = (R_0=6.0, D_lat=1.3, D_inf=2.0, α=0.5, D_imm=10.5, ρ=0.7)

# SEIT4L
println("Running PMMH for SEIT4L with $n_particles particles, $n_samples samples, thinning=$thinning...")
println("This will take a while (~1-2 hours)...")
model_seit4l = pmmh_seit4l(flu_tdc.obs, n_particles)
chain_full = sample(model_seit4l, MH(proposal...), n_samples; init_params, progress=true)

# Thin the chain
chain_pmcmc = chain_full[1:thinning:end]
println("Thinned from $(length(chain_full)) to $(length(chain_pmcmc)) samples")

output_path = datadir("pmcmc_seit4l_chain_precomputed.jld2")
println("Saving SEIT4L chain to $output_path")
@save output_path chain_pmcmc

println("\nSEIT4L chain summary:")
println(chain_pmcmc)

# SEITL
println("\nRunning PMMH for SEITL with $n_particles particles, $n_samples samples, thinning=$thinning...")
println("This will take a while (~1-2 hours)...")
model_seitl = pmmh_seitl(flu_tdc.obs, n_particles)
chain_seitl_full = sample(model_seitl, MH(proposal...), n_samples; init_params, progress=true)

# Thin the chain
chain_seitl = chain_seitl_full[1:thinning:end]
println("Thinned from $(length(chain_seitl_full)) to $(length(chain_seitl)) samples")

output_path_seitl = datadir("pmcmc_seitl_chain_precomputed.jld2")
println("Saving SEITL chain to $output_path_seitl")
@save output_path_seitl chain_seitl

println("\nSEITL chain summary:")
println(chain_seitl)

println("\nAll done!")
