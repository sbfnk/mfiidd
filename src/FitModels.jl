using DifferentialEquations
using Distributions
using DataFrames
using Plots
using Turing
using StatsPlots
using MCMCChains

"""
    FitModel

A structure that defines a model for fitting, compatible with Turing.jl patterns.
Contains:
- name: descriptive name of the model
- state_names: names of state variables
- theta_names: names of parameters
- priors: Dictionary of parameter names to Distribution objects
- simulate: function to simulate the model
- likelihood: function to calculate likelihood using Distribution objects
"""
struct FitModel
    name::String
    state_names::Vector{String}
    theta_names::Vector{String}
    priors::Dict{Symbol, Distribution}
    simulate::Function
    likelihood::Function
end

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
    sir_likelihood(data, trajectory)

Calculate likelihood of data given model trajectory using Distributions.jl.
Assumes observations follow Poisson distribution around true prevalence.
"""
function sir_likelihood(data, trajectory)
    # Initialize log-likelihood
    total_loglik = 0.0
    
    # Calculate likelihood for each observation
    for i in 1:nrow(data)
        # Find corresponding model point
        time_idx = findfirst(trajectory.time .≈ data.time[i])
        if time_idx !== nothing
            # Poisson likelihood: obs ~ Poisson(I)
            obs_dist = Poisson(max(trajectory.I[time_idx], 1e-6))  # Avoid zero rates
            total_loglik += logpdf(obs_dist, data.obs[i])
        end
    end
    
    return total_loglik
end

"""
    create_sir_deterministic()

Create deterministic SIR fitmodel object with Turing.jl compatible structure.
"""
function create_sir_deterministic()
    sir_name = "SIR with constant population size"
    sir_state_names = ["S", "I", "R"]
    sir_theta_names = ["R_0", "D_inf"]
    
    # Define priors using Distribution objects
    sir_priors = Dict(
        :R_0 => Uniform(1.0, 100.0),      # R_0 ~ Uniform(1, 100)
        :D_inf => Uniform(0.1, 30.0)      # D_inf ~ Uniform(0.1, 30)
    )
    
    return FitModel(
        sir_name,
        sir_state_names,
        sir_theta_names,
        sir_priors,
        sir_simulate_deterministic,
        sir_likelihood
    )
end

"""
    log_prior(fitmodel, theta)

Calculate log-prior density using Distribution objects.
"""
function log_prior(fitmodel, theta)
    total_logprior = 0.0
    for (param, prior_dist) in fitmodel.priors
        if haskey(theta, param)
            total_logprior += logpdf(prior_dist, theta[param])
        elseif hasproperty(theta, param)
            total_logprior += logpdf(prior_dist, getproperty(theta, param))
        else
            error("Parameter $param not found in theta")
        end
    end
    return total_logprior
end

"""
    log_likelihood(fitmodel, theta, init_state, data)

Calculate log-likelihood of trajectory with respect to data.
"""
function log_likelihood(fitmodel, theta, init_state, data)
    # Simulate the model
    trajectory = fitmodel.simulate(theta, init_state, data.time)
    
    # Calculate likelihood using the model's likelihood function
    return fitmodel.likelihood(data, trajectory)
end

"""
    log_posterior(fitmodel, theta, init_state, data)

Calculate log-posterior density for given parameters and initial state.
This follows Turing.jl patterns: log_posterior = log_prior + log_likelihood
"""
function log_posterior(fitmodel, theta, init_state, data)
    # Calculate log-prior
    logprior = log_prior(fitmodel, theta)
    
    # Calculate log-likelihood  
    loglik = log_likelihood(fitmodel, theta, init_state, data)
    
    # Return log-posterior
    return logprior + loglik
end

"""
    @model sir_model(data, init_state)

Turing.jl model definition for SIR model.
This demonstrates how to use the FitModel structure with Turing.jl.
"""
@model function sir_model(data, init_state)
    # Priors - using the same distributions as in our FitModel
    R_0 ~ Uniform(1.0, 100.0)
    D_inf ~ Uniform(0.1, 30.0)
    
    # Package parameters
    theta = (R_0=R_0, D_inf=D_inf)
    
    # Simulate model
    trajectory = sir_simulate_deterministic(theta, init_state, data.time)
    
    # Likelihood - observations follow Poisson distribution
    for i in 1:nrow(data)
        # Find corresponding time point in trajectory
        time_idx = findfirst(trajectory.time .≈ data.time[i])
        if time_idx !== nothing
            expected_cases = max(trajectory.I[time_idx], 1e-6)  # Avoid zero rates
            data.obs[i] ~ Poisson(expected_cases)
        end
    end
end

"""
    sample_posterior(fitmodel, data, init_state; n_samples=1000, sampler=NUTS())

Sample from posterior using Turing.jl.
"""
function sample_posterior(fitmodel, data, init_state; n_samples=1000, sampler=NUTS())
    # Create the Turing model
    model = sir_model(data, init_state)
    
    # Sample from posterior
    chain = sample(model, sampler, n_samples)
    
    return chain
end

"""
    r_traj_obs(fitmodel, theta, init_state, times)

Generate trajectory with simulated observations using Distribution objects.
"""
function r_traj_obs(fitmodel, theta, init_state, times)
    # Simulate the model
    trajectory = fitmodel.simulate(theta, init_state, times)
    
    # Generate observations at each time point using Poisson distribution
    obs = zeros(Int, length(times))
    for i in 1:length(times)
        expected_cases = max(trajectory.I[i], 1e-6)  # Avoid zero rates
        obs_dist = Poisson(expected_cases)
        obs[i] = rand(obs_dist)
    end
    
    # Add observations to trajectory
    trajectory.obs = obs
    
    return trajectory
end

"""
    plot_trajectory(trajectory; data=nothing, title="Model Trajectory")

Plot model trajectory, optionally with data overlay.
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

"""
    plot_fit(fitmodel, theta, init_state, data; n_replicates=1)

Plot model fit against data with optional replicates.
"""
function plot_fit(fitmodel, theta, init_state, data; n_replicates=1)
    p = plot(title="Model Fit")
    
    # Plot replicates
    for i in 1:n_replicates
        obs_traj = r_traj_obs(fitmodel, theta, init_state, data.time)
        if i == 1
            plot!(p, obs_traj.time, obs_traj.obs, alpha=0.6, color=:blue, label="Simulated", linewidth=1)
        else
            plot!(p, obs_traj.time, obs_traj.obs, alpha=0.6, color=:blue, label="", linewidth=1)
        end
    end
    
    # Plot data
    scatter!(p, data.time, data.obs, label="Observed", markersize=4, color=:black)
    
    xlabel!(p, "Time")
    ylabel!(p, "Observations")
    
    return p
end

"""
    plot_posterior_predictive(chain, data, init_state; n_samples=100)

Plot posterior predictive checks using samples from MCMC chain.
"""
function plot_posterior_predictive(chain, data, init_state; n_samples=100)
    p = plot(title="Posterior Predictive Check")
    
    # Extract parameter samples
    n_total = length(chain)
    sample_indices = rand(1:n_total, min(n_samples, n_total))
    
    for i in sample_indices
        # Extract parameters from chain
        R_0_sample = chain[:R_0][i]
        D_inf_sample = chain[:D_inf][i]
        theta_sample = Dict(:R_0 => R_0_sample, :D_inf => D_inf_sample)
        
        # Generate trajectory
        trajectory = sir_simulate_deterministic(theta_sample, init_state, data.time)
        
        # Plot with transparency
        if i == sample_indices[1]
            plot!(p, trajectory.time, trajectory.I, alpha=0.1, color=:blue, 
                  linewidth=0.5, label="Posterior samples")
        else
            plot!(p, trajectory.time, trajectory.I, alpha=0.1, color=:blue, 
                  linewidth=0.5, label="")
        end
    end
    
    # Plot data
    scatter!(p, data.time, data.obs, label="Observed", markersize=4, color=:black)
    
    xlabel!(p, "Time")
    ylabel!(p, "Cases")
    
    return p
end