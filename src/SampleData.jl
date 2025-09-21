using DataFrames
using Random

"""
    create_epi1_data()

Create sample epidemic data similar to the R package's epi1 dataset.
Generated from SIR model with R_0≈2.5, D_inf=2, N=1000, I0=1.
"""
function create_epi1_data()
    # Time points (weekly observations)
    times = 0:1:40
    
    # Simulated observations based on the original epi1 data pattern
    # These approximate the pattern from the R epi1 dataset
    obs = [1, 1, 1, 3, 4, 8, 12, 18, 25, 30, 35, 38, 40, 39, 35, 31, 
           26, 21, 17, 13, 11, 8, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1, 
           0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    return DataFrame(time = times, obs = obs)
end

"""
    create_epi2_data()

Create sample epidemic data with 10% reporting rate.
"""
function create_epi2_data()
    times = 0:1:40
    
    # Simulated observations with lower reporting (roughly 10% of epi1)
    obs = [0, 0, 0, 0, 1, 1, 1, 2, 3, 3, 4, 4, 4, 4, 4, 3, 
           3, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    return DataFrame(time = times, obs = obs)
end

"""
    load_models()

Create and return standard fitmodel objects.
"""
function load_models()
    sir_deter = create_sir_deterministic()
    
    return Dict(
        :sir_deter => sir_deter
    )
end

"""
    load_epi_data()

Load epidemic datasets.
"""
function load_epi_data()
    return Dict(
        :epi1 => create_epi1_data(),
        :epi2 => create_epi2_data()
    )
end