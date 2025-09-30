#!/usr/bin/env julia

# Test script to verify the Julia translation works with Turing.jl compatibility
println("Testing Turing.jl compatible Julia translation...")

try
    # Test loading modules
    include("src/FitModels.jl")
    include("src/SampleData.jl")
    println("✓ Successfully loaded modules")
    
    
    # Test simulation
    theta = (R_0=2.5, D_inf=2.0)
    init_state = [999.0, 1.0, 0.0]
    times = 0:1:20
    
    trajectory = sir_simulate_deterministic(theta, init_state, times)
    println("✓ Successfully simulated model trajectory")
    println("Trajectory size: ", size(trajectory))
    
    # Test data loading
    epi_data = load_epi_data()
    epi1 = epi_data[:epi1]
    println("✓ Successfully loaded epidemic data")
    println("Data size: ", size(epi1))
    
    # Test Turing.jl model and conditioning
    model = sir_model(epi1.obs, init_state, epi1.time)
    println("✓ Successfully created Turing model")
    
    # Test conditioning
    conditioned_model = model | theta  
    posterior_val = logpdf(conditioned_model, ())
    println("✓ Successfully calculated log posterior using conditioning: ", posterior_val)
    
    # Test that we can sample (just a few samples for testing)
    small_data_obs = epi1.obs[1:10]  # Use smaller dataset for faster testing  
    small_times = epi1.time[1:10]
    small_model = sir_model(small_data_obs, init_state, small_times)
    
    chain = sample(small_model, Prior(), 10)  # Sample from prior for quick test
    println("✓ Successfully sampled from Turing model")
    
    println("\n🎉 All tests passed! Turing.jl compatible Julia translation is working correctly.")
    
catch e
    println("❌ Error during testing: ", e)
    rethrow(e)
end