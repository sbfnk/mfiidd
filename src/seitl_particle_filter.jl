using Random
using Distributions
using StatsBase

"""
    gillespie_step_seitl!(rng, state, θ, dt)
    gillespie_step_seitl!(state, θ, dt)

Simulate SEITL for one time unit using the Gillespie algorithm.

The two-argument form draws from the global random number generator; pass an
`rng` explicitly when the caller needs to control randomness, as the particle
filter does.

# Arguments
- `rng`: Random number generator (defaults to the global one)
- `state`: Vector [S, E, I, T, L] (modified in place)
- `θ`: Parameter dictionary with keys :R_0, :D_lat, :D_inf, :α, :D_imm
- `dt`: Time step (typically 1.0 for daily)

# Returns
- `incidence`: Number of new cases (E→I transitions)
"""
function gillespie_step_seitl!(rng::AbstractRNG, state::Vector{Float64}, θ::Dict,
                               dt::Float64=1.0)
    β = θ[:R_0] / θ[:D_inf]
    ϵ = 1.0 / θ[:D_lat]
    ν = 1.0 / θ[:D_inf]
    τ = 1.0 / θ[:D_imm]
    α = θ[:α]

    # Stoichiometry: [S, E, I, T, L]
    stoich = [
        [-1, 1, 0, 0, 0],   # S → E (infection)
        [0, -1, 1, 0, 0],   # E → I (becoming infectious)
        [0, 0, -1, 1, 0],   # I → T (recovery)
        [1, 0, 0, -1, 0],   # T → S (immunity wanes)
        [0, 0, 0, -1, 1]    # T → L (long-term immunity)
    ]

    function rates(s)
        S, E, I, T, L = s
        N = S + E + I + T + L
        [β*S*I/N, ϵ*E, ν*I, (1-α)*τ*T, α*τ*T]
    end

    t, incidence = 0.0, 0
    while t < dt
        r = rates(state)
        total_rate = sum(r)
        total_rate ≤ 0 && break

        τ_wait = randexp(rng) / total_rate
        t + τ_wait > dt && break
        t += τ_wait

        # Select event
        cum, rnd, event = 0.0, rand(rng) * total_rate, 0
        for i in 1:5
            cum += r[i]
            if rnd ≤ cum
                event = i
                break
            end
        end

        # Apply transition
        for j in 1:5
            state[j] += stoich[event][j]
        end

        # E→I transitions count as new cases
        event == 2 && (incidence += 1)
    end

    return incidence
end

function gillespie_step_seitl!(state::Vector{Float64}, θ::Dict, dt::Float64=1.0)
    return gillespie_step_seitl!(Random.default_rng(), state, θ, dt)
end
