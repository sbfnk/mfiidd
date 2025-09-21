#!/usr/bin/env julia

# Test script to verify the Julia translation works with Turing.jl compatibility
println("Testing Turing.jl compatible Julia translation...")

try
    # Test loading modules
    include("src/FitModels.jl")
    include("src/SampleData.jl")
    println("✓ Successfully loaded modules")
    
    # Test creating models
    models = load_models()
    sir_deter = models[:sir_deter]
    println("✓ Successfully created SIR model")
    
    # Test model properties
    println("Model name: ", sir_deter.name)
    println("State names: ", sir_deter.state_names)
    println("Parameter names: ", sir_deter.theta_names)
    println("Priors: ", sir_deter.priors)
    
    # Test simulation
    theta = Dict(:R_0 => 2.5, :D_inf => 2.0)
    init_state = [999.0, 1.0, 0.0]
    times = 0:1:20
    
    trajectory = sir_deter.simulate(theta, init_state, times)
    println("✓ Successfully simulated model trajectory")
    println("Trajectory size: ", size(trajectory))
    
    # Test prior using new Distribution objects
    prior_val = log_prior(sir_deter, theta)
    println("✓ Successfully calculated log prior: ", prior_val)
    
    # Test data loading
    epi_data = load_epi_data()
    epi1 = epi_data[:epi1]
    println("✓ Successfully loaded epidemic data")
    println("Data size: ", size(epi1))
    
    # Test likelihood using new structure
    likelihood_val = log_likelihood(sir_deter, theta, init_state, epi1)
    println("✓ Successfully calculated log likelihood: ", likelihood_val)
    
    # Test posterior using new structure
    posterior_val = log_posterior(sir_deter, theta, init_state, epi1)
    println("✓ Successfully calculated log posterior: ", posterior_val)
    
    # Test observation generation
    obs_trajectory = r_traj_obs(sir_deter, theta, init_state, epi1.time)
    println("✓ Successfully generated observation trajectory")
    
    # Test Turing.jl integration
    println("Testing Turing.jl integration...")
    small_data = epi1[1:10, :]  # Use smaller dataset for faster testing
    
    # Test model compilation (don't run full sampling in test)
    model = sir_model(small_data, init_state)
    println("✓ Successfully created Turing model")
    
    # Test that we can sample (just a few samples for testing)
    chain = sample(model, Prior(), 10)  # Sample from prior for quick test
    println("✓ Successfully sampled from Turing model")
    
    println("\n🎉 All tests passed! Turing.jl compatible Julia translation is working correctly.")
    
catch e
    println("❌ Error during testing: ", e)
    rethrow(e)
end