#!/usr/bin/env julia

# Test script to verify the Julia translation works
println("Testing Julia translation of MFIIDD course materials...")

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
    
    # Test simulation
    theta = Dict(:R_0 => 2.5, :D_inf => 2.0)
    init_state = [999.0, 1.0, 0.0]
    times = 0:1:20
    
    trajectory = sir_deter.simulate(theta, init_state, times)
    println("✓ Successfully simulated model trajectory")
    println("Trajectory size: ", size(trajectory))
    
    # Test prior
    prior_val = sir_deter.d_prior(theta, log=true)
    println("✓ Successfully calculated prior: ", prior_val)
    
    # Test likelihood
    data_point = Dict(:obs => 18)
    model_point = Dict(:I => 31)
    likelihood_val = sir_deter.d_point_obs(data_point, model_point, theta, log=true)
    println("✓ Successfully calculated likelihood: ", likelihood_val)
    
    # Test data loading
    epi_data = load_epi_data()
    epi1 = epi_data[:epi1]
    println("✓ Successfully loaded epidemic data")
    println("Data size: ", size(epi1))
    
    # Test trajectory likelihood
    traj_likelihood = d_traj_obs(sir_deter, theta, init_state, epi1, log=true)
    println("✓ Successfully calculated trajectory likelihood: ", traj_likelihood)
    
    # Test posterior
    posterior_val = my_d_log_posterior(sir_deter, theta, init_state, epi1)
    println("✓ Successfully calculated posterior: ", posterior_val)
    
    # Test observation generation
    obs_trajectory = r_traj_obs(sir_deter, theta, init_state, epi1.time)
    println("✓ Successfully generated observation trajectory")
    
    println("\n🎉 All tests passed! Julia translation is working correctly.")
    
catch e
    println("❌ Error during testing: ", e)
    rethrow(e)
end