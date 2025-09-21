using DifferentialEquations
using Distributions
using DataFrames
using Plots

"""
    FitModel

A structure that defines a model for fitting, containing:
- name: descriptive name of the model
- state_names: names of state variables
- theta_names: names of parameters
- simulate: function to simulate the model
- d_prior: function to calculate prior density
- d_point_obs: function to calculate likelihood of a data point
- r_point_obs: function to generate observations from model
"""
struct FitModel
    name::String
    state_names::Vector{String}
    theta_names::Vector{String}
    simulate::Function
    d_prior::Function
    d_point_obs::Function
    r_point_obs::Function
end

"""
    sir_simulate_deterministic(theta, init_state, times)

Simulate deterministic SIR model with parameters theta, initial state init_state
at times specified in times vector.
"""
function sir_simulate_deterministic(theta, init_state, times)
    # Extract parameters
    R_0 = theta[:R_0]
    D_inf = theta[:D_inf]
    
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
    sir_prior(theta; log=false)

Calculate prior density for SIR model parameters.
Uses uniform priors: R_0 ~ U(1,100), D_inf ~ U(0,30)
"""
function sir_prior(theta; log=false)
    # Uniform prior on R_0: U[1,100]
    log_prior_R0 = logpdf(Uniform(1, 100), theta[:R_0])
    # Uniform prior on infectious period: U[0,30]
    log_prior_D_inf = logpdf(Uniform(0, 30), theta[:D_inf])
    
    log_sum = log_prior_R0 + log_prior_D_inf
    
    return log ? log_sum : exp(log_sum)
end

"""
    sir_point_likelihood(data_point, model_point, theta; log=false)

Calculate likelihood of a single data point given model state.
Assumes observations follow Poisson distribution around true prevalence.
"""
function sir_point_likelihood(data_point, model_point, theta; log=false)
    # The prevalence is observed through a Poisson process
    return logpdf(Poisson(model_point[:I]), data_point[:obs]) * (log ? 1 : exp(1))
end

"""
    sir_generate_obs_point(model_point, theta)

Generate random observation from model state.
"""
function sir_generate_obs_point(model_point, theta)
    # The prevalence is observed through a Poisson process
    obs_point = rand(Poisson(model_point[:I]))
    return Dict(:obs => obs_point)
end

"""
    create_sir_deterministic()

Create deterministic SIR fitmodel object.
"""
function create_sir_deterministic()
    sir_name = "SIR with constant population size"
    sir_state_names = ["S", "I", "R"]
    sir_theta_names = ["R_0", "D_inf"]
    
    return FitModel(
        sir_name,
        sir_state_names,
        sir_theta_names,
        sir_simulate_deterministic,
        sir_prior,
        sir_point_likelihood,
        sir_generate_obs_point
    )
end

"""
    d_traj_obs(fitmodel, theta, init_state, data; log=true)

Calculate log-likelihood of trajectory with respect to data.
"""
function d_traj_obs(fitmodel, theta, init_state, data; log=true)
    # 1. Simulate the model
    trajectory = fitmodel.simulate(theta, init_state, data.time)
    
    # 2. Calculate likelihood at each data point
    log_likelihood = 0.0
    for i in 1:nrow(data)
        data_point = Dict(:obs => data.obs[i])
        model_point = Dict(:I => trajectory.I[i])
        log_likelihood += fitmodel.d_point_obs(data_point, model_point, theta, log=true)
    end
    
    return log ? log_likelihood : exp(log_likelihood)
end

"""
    my_d_log_posterior(fitmodel, theta, init_state, data)

Calculate log-posterior density for given parameters and initial state.
"""
function my_d_log_posterior(fitmodel, theta, init_state, data)
    # Calculate log-prior
    log_prior = fitmodel.d_prior(theta, log=true)
    
    # Calculate log-likelihood
    log_likelihood = d_traj_obs(fitmodel, theta, init_state, data, log=true)
    
    # Calculate log-posterior
    log_posterior = log_prior + log_likelihood
    
    return log_posterior
end

"""
    r_traj_obs(fitmodel, theta, init_state, times)

Generate trajectory with simulated observations.
"""
function r_traj_obs(fitmodel, theta, init_state, times)
    # 1. Simulate the model
    trajectory = fitmodel.simulate(theta, init_state, times)
    
    # 2. Generate observations at each time point
    obs = zeros(Int, length(times))
    for i in 1:length(times)
        model_point = Dict(:I => trajectory.I[i])
        obs_point = fitmodel.r_point_obs(model_point, theta)
        obs[i] = obs_point[:obs]
    end
    
    # 3. Add observations to trajectory
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