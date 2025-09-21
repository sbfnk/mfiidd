#!/usr/bin/env julia

# Demonstration of elegant Turing.jl integration for Bayesian infectious disease modeling
println("🦠 Demonstrating modern Bayesian infectious disease modeling with Turing.jl")
println("=" ^ 80)

include("src/FitModels.jl")
include("src/SampleData.jl")

# Load the model and data
models = load_models()
sir_deter = models[:sir_deter]
epi_data = load_epi_data()
epi1 = epi_data[:epi1][1:25, :]  # Use subset for faster demo

println("📊 Model: ", sir_deter.name)
println("📊 Data points: ", nrow(epi1))
println()

# Show elegant prior specification using Distributions.jl
println("🎯 Prior distributions (using Distributions.jl):")
for (param, prior) in sir_deter.priors
    println("  $param ~ $prior")
end
println()

# Traditional approach: manual posterior exploration
println("📈 Traditional approach: Grid search over parameter space")
init_state = [999.0, 1.0, 0.0]

R_0_values = 1.5:0.5:4.0
posteriors = []
for R_0 in R_0_values
    theta = Dict(:R_0 => R_0, :D_inf => 2.0)
    posterior = log_posterior(sir_deter, theta, init_state, epi1)
    push!(posteriors, posterior)
end

best_idx = argmax(posteriors)
best_R_0 = R_0_values[best_idx]
println("  Best R_0 from grid search: ", best_R_0)
println()

# Modern approach: Turing.jl MCMC
println("🚀 Modern approach: Bayesian inference with Turing.jl")
println("  Sampling from posterior using NUTS sampler...")

# Sample from posterior
chain = sample_posterior(sir_deter, epi1, init_state, n_samples=500, sampler=NUTS())

println("✓ MCMC sampling complete!")
println()

# Extract results
R_0_samples = chain[:R_0]
D_inf_samples = chain[:D_inf]

println("📊 Posterior summary:")
println("  R_0: mean = ", round(mean(R_0_samples), digits=3), 
        ", std = ", round(std(R_0_samples), digits=3))
println("  D_inf: mean = ", round(mean(D_inf_samples), digits=3), 
        ", std = ", round(std(D_inf_samples), digits=3))
println()

# Compare approaches
println("🔍 Comparison:")
println("  Grid search best R_0: ", best_R_0)
println("  MCMC posterior mean R_0: ", round(mean(R_0_samples), digits=3))
println("  Difference: ", round(abs(best_R_0 - mean(R_0_samples)), digits=3))
println()

# Show the power of the Bayesian approach
println("🎉 Benefits of the Turing.jl approach:")
println("  ✓ Elegant prior specification using Distribution objects")
println("  ✓ Natural ~ syntax for likelihood specification")  
println("  ✓ Automatic MCMC sampling with modern algorithms (NUTS)")
println("  ✓ Built-in diagnostics and convergence checks")
println("  ✓ Posterior predictive checks")
println("  ✓ Easy to extend to more complex models")
println()

# Demonstrate posterior predictive checks
println("📈 Generating posterior predictive samples for model validation...")
using Plots
p = plot_posterior_predictive(chain, epi1, init_state, n_samples=20)
println("✓ Posterior predictive plot generated (see plot window)")
println()

# Demonstrate easy model modification
println("🔧 Easy model modification with Distribution objects:")
println("  Original R_0 prior: ", sir_deter.priors[:R_0])

# Create model with informative prior
function create_informative_sir()
    return FitModel(
        "SIR with informative priors",
        sir_deter.state_names,
        sir_deter.theta_names,
        Dict(:R_0 => Normal(2.5, 0.3), :D_inf => Gamma(4, 0.5)),
        sir_deter.simulate,
        sir_deter.likelihood
    )
end

sir_informative = create_informative_sir()
println("  Informative R_0 prior: ", sir_informative.priors[:R_0])

# Compare posteriors
theta_test = Dict(:R_0 => 2.5, :D_inf => 2.0)
posterior_flat = log_posterior(sir_deter, theta_test, init_state, epi1)
posterior_informative = log_posterior(sir_informative, theta_test, init_state, epi1)

println("  Log posterior (flat prior): ", round(posterior_flat, digits=2))
println("  Log posterior (informative prior): ", round(posterior_informative, digits=2))
println("  Prior influence: ", round(posterior_informative - posterior_flat, digits=2))
println()

println("🎯 Ready for advanced Bayesian infectious disease modeling!")
println("   Next steps: Extend to stochastic models, multi-strain dynamics, etc.")
println("=" ^ 80)