# Adaptive PMMH using empirical covariance proposals
# This implements Adaptive Metropolis (Haario et al. 2001) for PMMH
#
# Run with: julia --project=. scripts/adaptive_pmmh.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
using Distributions
using DataFrames
using Turing
using CSV
using DrWatson
using StatsBase
using LinearAlgebra
using MCMCChains

Random.seed!(42)

# Load data
flu_tdc = CSV.read(datadir("flu_tdc_1971.csv"), DataFrame)

# =============================================================================
# Particle filter for SEIT4L (same as before)
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
# Log-posterior function
# =============================================================================

# Prior bounds
const BOUNDS = (
    R_0 = (1.0, 50.0),
    D_lat = (0.1, 10.0),
    D_inf = (0.1, 15.0),
    α = (0.0, 1.0),
    D_imm = (1.0, 50.0),
    ρ = (0.0, 1.0)
)
const PARAM_NAMES = [:R_0, :D_lat, :D_inf, :α, :D_imm, :ρ]

function log_prior(x)
    # Uniform priors - check bounds
    for (i, p) in enumerate(PARAM_NAMES)
        lo, hi = BOUNDS[p]
        if x[i] < lo || x[i] > hi
            return -Inf
        end
    end
    # Sum of log(1/(hi-lo)) for uniform priors
    lp = 0.0
    for p in PARAM_NAMES
        lo, hi = BOUNDS[p]
        lp -= log(hi - lo)
    end
    return lp
end

function log_posterior(x, obs, n_particles)
    lp = log_prior(x)
    lp == -Inf && return -Inf

    θ = Dict(p => x[i] for (i, p) in enumerate(PARAM_NAMES))
    ll = particle_filter_seit4l(θ, obs, n_particles)
    return lp + ll
end

# =============================================================================
# Adaptive Metropolis (Haario et al. 2001)
# =============================================================================

function adaptive_mcmc(obs, n_particles, n_samples;
                       init_x=nothing,
                       init_Σ=nothing,
                       adapt_start=100,
                       ε=1e-6)
    d = 6  # number of parameters

    # Initial state
    if init_x === nothing
        init_x = [6.0, 1.3, 2.0, 0.5, 10.5, 0.7]  # near posterior mode
    end

    # Initial proposal covariance (diagonal, will adapt)
    if init_Σ === nothing
        init_Σ = diagm([4.0, 0.3, 1.5, 0.05, 2.0, 0.05].^2)
    end

    # Scaling factor (optimal for d-dimensional Gaussian: 2.38²/d)
    s_d = 2.38^2 / d

    # Storage
    samples = zeros(n_samples, d)
    log_posts = zeros(n_samples)

    # Running mean and covariance for adaptation
    x_mean = copy(init_x)
    x_cov = copy(init_Σ)

    # Current state
    x = copy(init_x)
    lp = log_posterior(x, obs, n_particles)

    n_accept = 0

    println("Starting Adaptive MCMC with $n_samples samples...")
    t_start = time()

    for i in 1:n_samples
        # Propose using adapted covariance (+ small diagonal for stability)
        Σ_prop = s_d * x_cov + ε * I(d)

        # Ensure positive definite
        Σ_prop = (Σ_prop + Σ_prop') / 2  # symmetrize

        x_prop = x + rand(MvNormal(zeros(d), Σ_prop))
        lp_prop = log_posterior(x_prop, obs, n_particles)

        # Accept/reject
        if log(rand()) < lp_prop - lp
            x = x_prop
            lp = lp_prop
            n_accept += 1
        end

        # Store sample
        samples[i, :] = x
        log_posts[i] = lp

        # Update running statistics for adaptation (after adapt_start)
        if i >= adapt_start
            # Welford's online algorithm for mean and covariance
            n = i - adapt_start + 1
            x_old_mean = copy(x_mean)
            x_mean = x_old_mean + (x - x_old_mean) / n
            if n > 1
                x_cov = ((n - 2) * x_cov + (x - x_old_mean) * (x - x_mean)') / (n - 1)
            end
        end

        # Progress
        if i % 1000 == 0
            elapsed = time() - t_start
            rate = i / elapsed
            accept_rate = 100 * n_accept / i
            println("  $i/$n_samples ($(round(rate, digits=1))/s, accept=$(round(accept_rate, digits=1))%)")
        end
    end

    elapsed = time() - t_start
    println("\nCompleted in $(round(elapsed/60, digits=1)) minutes")
    println("Overall acceptance rate: $(round(100*n_accept/n_samples, digits=1))%")

    return samples, log_posts, x_cov
end

# =============================================================================
# Run
# =============================================================================

n_particles = 256
n_samples = 20000  # Adaptive should need fewer samples

# Use empirical covariance from previous run as starting point
# (This gives adaptation a head start)
Σ_init = [
    13.52  0.5361  2.216  0.01045  1.378  0.0004056;
    0.5361  0.1486  -0.2839  -0.001504  0.1571  0.002234;
    2.216  -0.2839  2.132  0.01455  -0.2214  -0.01215;
    0.01045  -0.001504  0.01455  0.0003732  -0.005007  -0.0001918;
    1.378  0.1571  -0.2214  -0.005007  1.566  0.001945;
    0.0004056  0.002234  -0.01215  -0.0001918  0.001945  0.0003544
]

samples, log_posts, final_cov = adaptive_mcmc(
    flu_tdc.obs, n_particles, n_samples;
    init_Σ = Σ_init,
    adapt_start = 500  # Start adapting after 500 samples
)

# Convert to MCMCChains for diagnostics
chain = Chains(samples, PARAM_NAMES)
println("\n=== Chain Summary ===")
println(describe(chain)[1])

println("\n=== Effective Sample Size ===")
ess_result = ess(chain)
for p in PARAM_NAMES
    println("  $p: $(round(ess_result[p, :ess], digits=1))")
end

# Save
df = DataFrame(samples, PARAM_NAMES)
df.logjoint = log_posts
output_path = datadir("adaptive_pmcmc_seit4l_samples.csv")
CSV.write(output_path, df)
println("\nSaved $(nrow(df)) samples to $output_path")

println("\n=== Final adapted covariance (for future runs) ===")
println("Σ_adapted = [")
for i in 1:6
    vals = [round(final_cov[i,j], sigdigits=4) for j in 1:6]
    println("  ", join(vals, "  "), ";")
end
println("]")
