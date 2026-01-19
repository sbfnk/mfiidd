# Proper bootstrap particle filter with resampling
# Replicates the R SSM approach for calibration

using Random, Distributions, DataFrames, CSV, DrWatson, Plots, StatsBase

Random.seed!(42)
flu_tdc = CSV.read(datadir("flu_tdc_1971.csv"), DataFrame)

# SEIT4L transition rates
function seit4l_rates(state, θ)
    S, E, I, T1, T2, T3, T4, L = state
    R_0, D_lat, D_inf, α, D_imm = θ[:R_0], θ[:D_lat], θ[:D_inf], θ[:α], θ[:D_imm]

    β = R_0 / D_inf
    ϵ = 1.0 / D_lat
    ν = 1.0 / D_inf
    τ = 4.0 / D_imm
    N = S + E + I + T1 + T2 + T3 + T4 + L

    return [
        β * S * I / N,      # S → E
        ϵ * E,              # E → I
        ν * I,              # I → T1
        τ * T1,             # T1 → T2
        τ * T2,             # T2 → T3
        τ * T3,             # T3 → T4
        (1 - α) * τ * T4,   # T4 → S
        α * τ * T4          # T4 → L
    ]
end

# Stoichiometry for SEIT4L (without incidence tracking)
const STOICH = [
    -1  1  0  0  0  0  0  0;   # S → E
     0 -1  1  0  0  0  0  0;   # E → I
     0  0 -1  1  0  0  0  0;   # I → T1
     0  0  0 -1  1  0  0  0;   # T1 → T2
     0  0  0  0 -1  1  0  0;   # T2 → T3
     0  0  0  0  0 -1  1  0;   # T3 → T4
     1  0  0  0  0  0 -1  0;   # T4 → S
     0  0  0  0  0  0 -1  1    # T4 → L
]

# Simulate one step of Gillespie algorithm from t to t+dt, return new state and incidence
function gillespie_step!(state, θ, dt)
    t = 0.0
    incidence = 0

    while t < dt
        rates = seit4l_rates(state, θ)
        total_rate = sum(rates)

        if total_rate ≤ 0
            break
        end

        # Time to next event
        τ = randexp() / total_rate
        if t + τ > dt
            break
        end
        t += τ

        # Which event?
        r = rand() * total_rate
        cumsum = 0.0
        event = 0
        for i in 1:length(rates)
            cumsum += rates[i]
            if r ≤ cumsum
                event = i
                break
            end
        end

        # Apply stoichiometry
        for j in 1:8
            state[j] += STOICH[event, j]
        end

        # Track incidence (E → I transition, event 2)
        if event == 2
            incidence += 1
        end
    end

    return incidence
end

"""
Bootstrap particle filter with resampling.
Returns log-likelihood estimate and whether particle depletion occurred.
"""
function particle_filter(θ, obs, n_particles; init_state=[279.0, 0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0])
    n_obs = length(obs)
    ρ = θ[:ρ]

    # Initialize particles (each is a copy of initial state)
    particles = [copy(init_state) for _ in 1:n_particles]

    log_lik = 0.0
    depleted = false

    for t in 1:n_obs
        # Propagate particles forward by 1 time unit
        incidences = zeros(Int, n_particles)
        for i in 1:n_particles
            incidences[i] = gillespie_step!(particles[i], θ, 1.0)
        end

        # Compute weights (likelihood of observation given each particle)
        log_weights = zeros(n_particles)
        for i in 1:n_particles
            λ = max(ρ * incidences[i], 1e-10)
            log_weights[i] = logpdf(Poisson(λ), obs[t])
        end

        # Normalize weights (log-sum-exp trick)
        max_lw = maximum(log_weights)
        weights = exp.(log_weights .- max_lw)
        sum_weights = sum(weights)

        # Update log-likelihood
        log_lik += max_lw + log(sum_weights / n_particles)

        # Normalize for resampling
        weights ./= sum_weights

        # Check for particle depletion (effective sample size)
        ess = 1.0 / sum(weights.^2)
        if ess < n_particles / 2
            # Resample
            indices = wsample(1:n_particles, weights, n_particles)
            particles = [copy(particles[i]) for i in indices]

            # Check if all particles collapsed to one
            if length(unique(indices)) == 1
                depleted = true
            end
        end
    end

    return log_lik, depleted
end

"""
Run calibration: compute mean, SD of log-likelihood and depletion rate for different particle counts.
"""
function calibrate(θ, obs, particle_counts; n_reps=100)
    results = DataFrame(
        n_particles = Int[],
        mean_ll = Float64[],
        sd_ll = Float64[],
        prop_depleted = Float64[],
        time_per_iter = Float64[]
    )

    for n in particle_counts
        println("Testing n_particles = $n...")

        logliks = Float64[]
        n_depleted = 0

        t_start = time()
        for _ in 1:n_reps
            ll, depleted = particle_filter(θ, obs, n)
            push!(logliks, ll)
            if depleted
                n_depleted += 1
            end
        end
        t_elapsed = time() - t_start

        push!(results, (
            n_particles = n,
            mean_ll = mean(logliks),
            sd_ll = std(logliks),
            prop_depleted = n_depleted / n_reps,
            time_per_iter = t_elapsed / n_reps
        ))

        println("  mean=$(round(mean(logliks), digits=1)), SD=$(round(std(logliks), digits=2)), depleted=$(round(100*n_depleted/n_reps, digits=1))%")
    end

    return results
end

# Reference parameters (from deterministic fit)
θ_ref = Dict(
    :R_0 => 6.0,
    :D_lat => 1.3,
    :D_inf => 2.0,
    :α => 0.5,
    :D_imm => 10.5,
    :ρ => 0.7
)

# Run calibration
particle_counts = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
results = calibrate(θ_ref, flu_tdc.obs, particle_counts, n_reps=100)

println("\nResults:")
println(results)

# Save results
CSV.write(datadir("particle_calibration.csv"), results)
println("\nSaved to data/particle_calibration.csv")

# Create 4-panel plot like the R version
p1 = plot(results.n_particles, results.mean_ll,
    xlabel="Number of particles", ylabel="Mean log-likelihood",
    title="mean", legend=false, marker=:circle, linewidth=2)
vline!(p1, [128], color=:red, linestyle=:dash)

p2 = plot(results.n_particles, results.prop_depleted,
    xlabel="Number of particles", ylabel="Proportion",
    title="prop. of samples with particle depletion", legend=false, marker=:circle, linewidth=2)
vline!(p2, [128], color=:red, linestyle=:dash)

p3 = plot(results.n_particles, results.sd_ll,
    xlabel="Number of particles", ylabel="SD",
    title="sd", legend=false, marker=:circle, linewidth=2)
hline!(p3, [1, 3], color=:grey, linestyle=:dot, alpha=0.5)
vline!(p3, [128], color=:red, linestyle=:dash)

p4 = plot(results.n_particles, results.time_per_iter .* 10000,
    xlabel="Number of particles", ylabel="Time (seconds)",
    title="time 10000 iter", legend=false, marker=:circle, linewidth=2)
vline!(p4, [128], color=:red, linestyle=:dash)

p = plot(p1, p2, p3, p4, layout=(2,2), size=(800, 600))
savefig(p, datadir("particle_calibration.png"))
println("Saved plot to data/particle_calibration.png")
