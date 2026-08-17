# Shared setup for the PMMH chain generation scripts.
#
# The priors below must match those in sessions/pmcmc.qmd. The session presents
# the saved chains as the posteriors of the model it defines, so a prior changed
# in one place and not the other silently puts the wrong posterior on the page.

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
using ForwardDiff

using MFIIDD

# Sampler settings shared by every chain we generate
const N_PARTICLES = 256
const N_WARMUP = 50_000     # RAM adapts here, and these draws are discarded
const N_SAMPLES = 450_000   # kept, then thinned
const THINNING = 50         # → 9000 final samples

const PARAMETERS = [:R_0, :D_lat, :D_inf, :α, :D_imm, :ρ]

"""
    pmmh(obs, n_particles, particle_filter)

PMMH model: weakly informative priors on the six parameters, with the
log-likelihood estimated by `particle_filter` and added through `@addlogprob!`.

`ForwardDiff.value` strips the Duals that Turing's gradient-based initialisation
probe pushes through the model. The particle filter is not differentiable, and
its resampling step errors on Dual-valued weights.
"""
@model function pmmh(obs, n_particles, particle_filter)
    R_0 ~ truncated(Normal(3.0, 2.0), lower=1.0)
    D_lat ~ truncated(Normal(2.0, 1.0), lower=0.5)
    D_inf ~ truncated(Normal(3.0, 2.0), lower=0.5)
    α ~ Beta(2, 2)
    D_imm ~ truncated(Normal(15.0, 10.0), lower=1.0)
    ρ ~ Beta(2, 2)

    θ = Dict(:R_0 => ForwardDiff.value(R_0),
             :D_lat => ForwardDiff.value(D_lat),
             :D_inf => ForwardDiff.value(D_inf),
             :α => ForwardDiff.value(α),
             :D_imm => ForwardDiff.value(D_imm),
             :ρ => ForwardDiff.value(ρ))

    log_lik = particle_filter(θ, obs, n_particles)
    Turing.@addlogprob! log_lik
end

pmmh_seitl(obs, n_particles) = pmmh(obs, n_particles, particle_filter_seitl)
pmmh_seit4l(obs, n_particles) = pmmh(obs, n_particles, particle_filter_seit4l)

"""
    flu_observations()

Daily incidence from the 1971 Tristan da Cunha outbreak.
"""
flu_observations() = CSV.read(datadir("flu_tdc_1971.csv"), DataFrame).obs

"""
    chain_frame(chain)

Return `chain` as a data frame with one row per iteration and one column per
parameter, dropping the iteration and chain indices. Going through
`DataFrame(chain)` keeps this working across chain types rather than reaching
into internal fields whose layout is not part of the API.
"""
function chain_frame(chain)
    df = DataFrame(chain)
    return select(df, Not(intersect(["iteration", "iter", "chain"], names(df))))
end

"""
    save_chain_csv(chain, path)

Write the chain to `path`, one row per retained iteration.
"""
save_chain_csv(chain, path) = CSV.write(path, chain_frame(chain))

"""
    acceptance_rate(chain)

Proportion of iterations at which the sampler moved. A Metropolis chain repeats
the previous parameter vector whenever a proposal is rejected, so the iterations
where the parameters changed are the accepted ones.
"""
function acceptance_rate(chain)
    x = chain_frame(chain)[!, first(PARAMETERS)]
    return mean(x[2:end] .!= x[1:end-1])
end

"""
    run_pmmh(model, name; n_warmup, n_samples, thinning)

Sample `model` with Robust Adaptive Metropolis, reporting the acceptance rate of
the kept draws, then thin. Returns the thinned chain.

`num_warmup` is what makes RAM adaptive: the sampler only updates its proposal
covariance in the warmup phase, and those draws are discarded rather than kept.
Without it the proposal stays at the identity matrix and acceptance sits near
0.5% instead of the 23% RAM aims for. `adapting_externalsampler` is what carries
the warmup phase through to the sampler; see `src/adapting_sampler.jl`.
"""
function run_pmmh(model, name; n_warmup=N_WARMUP, n_samples=N_SAMPLES,
                  thinning=THINNING)
    println("=" ^ 60)
    println("Running PMMH for $name with RAM")
    println("  Particles: $N_PARTICLES")
    println("  Warmup (adaptation, discarded): $n_warmup")
    println("  Samples kept: $n_samples")
    println("  Thinning: $thinning")
    println("  Final samples: $(n_samples ÷ thinning)")
    println("=" ^ 60)

    t_start = time()
    chain_full = sample(
        model,
        adapting_externalsampler(AdvancedMH.RobustAdaptiveMetropolis()),
        n_samples;
        num_warmup=n_warmup,
        check_model=false,
        progress=true
    )
    t_elapsed = time() - t_start
    println("\n$name sampling took $(round(t_elapsed/60, digits=1)) minutes")
    println("Acceptance rate: $(round(acceptance_rate(chain_full) * 100, digits=1))%")

    chain = chain_full[1:thinning:end]
    println("After thinning: $(size(chain, 1)) samples")
    return chain
end

"""
    print_diagnostics(chain, name)

Print the posterior summary and credible intervals, from the same six columns
the session reads back out of the saved CSV.
"""
function print_diagnostics(chain, name)
    df = chain_frame(chain)
    mcmc_chain = Chains(Matrix(df[:, PARAMETERS]), PARAMETERS)

    println("\n$name summary statistics:")
    show(stdout, MIME("text/plain"), summarystats(mcmc_chain))
    println("\n\n$name quantiles:")
    show(stdout, MIME("text/plain"), quantile(mcmc_chain))
    println()
end
