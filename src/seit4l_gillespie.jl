using Random

"""
    gillespie_step(rng, state, θ, dt=1.0)

Simulate SEIT4L for `dt` time units using the Gillespie algorithm.

# Arguments
- `rng`: Random number generator
- `state`: Vector [S, E, I, T1, T2, T3, T4, L]
- `θ`: Parameter dictionary with keys :R_0, :D_lat, :D_inf, :α, :D_imm
- `dt`: Length of the simulation interval (default 1.0, i.e. one day)

# Returns
- `(new_state, incidence)`: Updated state and number of new cases over `dt`
"""
function gillespie_step(
    rng::AbstractRNG,
    state::Vector{Float64},
    θ::Dict,
    dt::Float64 = 1.0,
)
    β = θ[:R_0] / θ[:D_inf]
    ϵ = 1.0 / θ[:D_lat]
    ν = 1.0 / θ[:D_inf]
    τ = 4.0 / θ[:D_imm]
    α = θ[:α]

    # Copy state for modification
    s = copy(state)

    # Stoichiometry: how each transition changes [S, E, I, T1, T2, T3, T4, L]
    stoich = [
        [-1, 1, 0, 0, 0, 0, 0, 0],   # S → E (infection)
        [0, -1, 1, 0, 0, 0, 0, 0],   # E → I (becoming infectious)
        [0, 0, -1, 1, 0, 0, 0, 0],   # I → T1 (recovery)
        [0, 0, 0, -1, 1, 0, 0, 0],   # T1 → T2
        [0, 0, 0, 0, -1, 1, 0, 0],   # T2 → T3
        [0, 0, 0, 0, 0, -1, 1, 0],   # T3 → T4
        [1, 0, 0, 0, 0, 0, -1, 0],   # T4 → S (immunity wanes)
        [0, 0, 0, 0, 0, 0, -1, 1],    # T4 → L (long-term immunity)
    ]

    function rates(s)
        S, E, I, T1, T2, T3, T4, L = s
        N = S + E + I + T1 + T2 + T3 + T4 + L
        [β*S*I/N, ϵ*E, ν*I, τ*T1, τ*T2, τ*T3, (1-α)*τ*T4, α*τ*T4]
    end

    # Simulate up to `dt` time units
    t, daily_inc = 0.0, 0
    while t < dt
        r = rates(s)
        total_rate = sum(r)
        total_rate ≤ 0 && break

        # Time to next event
        τ_wait = randexp(rng) / total_rate
        t + τ_wait > dt && break
        t += τ_wait

        # Select which event occurs
        cum, rnd, event = 0.0, rand(rng) * total_rate, 0
        for i in 1:8
            cum += r[i]
            if rnd ≤ cum
                event = i
                break
            end
        end

        # Apply the transition
        for j in 1:8
            s[j] += stoich[event][j]
        end

        # E → I transitions count as new cases
        event == 2 && (daily_inc += 1)
    end

    return s, daily_inc
end

"""
    gillespie_step_seit4l!(state, θ, dt)

In-place version for simple bootstrap filter (no RNG argument, uses global RNG).
"""
function gillespie_step_seit4l!(state::Vector{Float64}, θ::Dict, dt::Float64 = 1.0)
    new_state, inc = gillespie_step(Random.default_rng(), state, θ, dt)
    for i in 1:8
        state[i] = new_state[i]
    end
    return inc
end
