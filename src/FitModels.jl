using DifferentialEquations
using ModellingToolkit
using Distributions
using DataFrames
using Plots
using Turing
using StatsPlots
using MCMCChains

# Core SIR simulation and Turing.jl model
# Use the sir_model for all Bayesian inference - it's the single source of truth

"""
    sir_simulate_deterministic(theta, init_state, times)

Simulate deterministic SIR model with parameters theta, initial state init_state
at times specified in times vector.
"""
function sir_simulate_deterministic(theta, init_state, times)
    # Extract parameters - support both Dict and NamedTuple
    R_0 = haskey(theta, :R_0) ? theta[:R_0] : theta.R_0
    D_inf = haskey(theta, :D_inf) ? theta[:D_inf] : theta.D_inf
    
    # Define ODE function
    function sir_ode!(du, u, p, t)
        S, I, R = u
        beta, nu = p
        N = S + I + R
        
        du[1] = -beta * S * I / N  # dS/dt
        du[2] = beta * S * I / N - nu * I  # dI/dt
        du[3] = nu * I  # dR/dt
    end
    
    # Parameters for ODE
    beta = R_0 / D_inf
    nu = 1 / D_inf
    params = [beta, nu]
    
    # Solve ODE
    prob = ODEProblem(sir_ode!, init_state, (times[1], times[end]), params)
    sol = solve(prob, Tsit5(), saveat=times)
    
    # Convert to DataFrame
    df = DataFrame(
        time = times,
        S = [sol(t)[1] for t in times],
        I = [sol(t)[2] for t in times], 
        R = [sol(t)[3] for t in times]
    )
    
    return df
end



"""
    @model sir_model(data, init_state, times)

Turing.jl model definition for SIR model - the single source of truth for Bayesian inference.
Use conditioning to evaluate at specific parameters: model | (R_0 = 2.5, D_inf = 2.0)
"""
@model function sir_model(data, init_state, times)
    # Priors
    R_0 ~ Uniform(1.0, 100.0)
    D_inf ~ Uniform(0.1, 30.0)
    
    # SIR ODE simulation inline
    function sir_ode!(du, u, p, t)
        S, I, R = u
        beta, nu = p
        N = S + I + R
        du[1] = -beta * S * I / N  # dS/dt
        du[2] = beta * S * I / N - nu * I  # dI/dt
        du[3] = nu * I  # dR/dt
    end
    
    # Solve ODE
    beta = R_0 / D_inf
    nu = 1 / D_inf
    prob = ODEProblem(sir_ode!, init_state, (times[1], times[end]), [beta, nu])
    sol = solve(prob, Tsit5(), saveat=times)
    
    # Likelihood - observations follow Poisson distribution
    for i in 1:length(data)
        expected_cases = max(sol.u[i][2], 1e-6)  # I compartment, avoid zero rates
        data[i] ~ Poisson(expected_cases)
    end
end


"""
    plot_trajectory(trajectory; data=nothing, title="Model Trajectory")

Plot SIR model trajectory, optionally with data overlay.
"""
function plot_trajectory(trajectory; data=nothing, title="Model Trajectory")
    p = plot(trajectory.time, trajectory.S, label="S", linewidth=2, title=title)
    plot!(p, trajectory.time, trajectory.I, label="I", linewidth=2)
    plot!(p, trajectory.time, trajectory.R, label="R", linewidth=2)
    
    if data !== nothing
        scatter!(p, data.time, data.obs, label="Observed", markersize=4, color=:black)
    end
    
    xlabel!(p, "Time")
    ylabel!(p, "Population")
    
    return p
end

function sir_ode!(du, u, p, t)
    S, I, R = u
    R_0, D_inf = p
    N = S + I + R

    β = R_0 / D_inf
    ν = 1 / D_inf

    du[1] = -β * S * I / N  # dS/dt
    du[2] = β * S * I / N - ν * I  # dI/dt
    du[3] = ν * I  # dR/dt
end
