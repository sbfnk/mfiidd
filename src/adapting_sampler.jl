using Turing: Turing, externalsampler
using Turing.DynamicPPL: DynamicPPL
using Turing.AbstractMCMC: AbstractMCMC
using Random: AbstractRNG

"""
    AdaptingExternalSampler(inner)

Wrapper that lets an adaptive external sampler actually adapt.

Turing's `externalsampler` forwards `AbstractMCMC.step` to the sampler it wraps
but not `AbstractMCMC.step_warmup`. AbstractMCMC then falls back to `step` for
warmup, and for an adaptive sampler that is the branch which does no adaptation.
`AdvancedMH.RobustAdaptiveMetropolis` only updates its proposal covariance from
`step_warmup`, so through a plain `externalsampler` it runs as a fixed random
walk with whatever proposal it started from, however long you sample. On the
SEITL model that is the difference between 0.5% and 23% acceptance.

Wrapping in our own sampler type lets us define the missing forwarding on a type
we own. Use [`adapting_externalsampler`](@ref) to construct one, and pass
`num_warmup` to `sample`; without warmup iterations there is nothing to adapt in.
"""
struct AdaptingExternalSampler{S<:Turing.Inference.ExternalSampler} <:
       AbstractMCMC.AbstractSampler
    inner::S
end

"""
    adapting_externalsampler(sampler; kwargs...)

Wrap `sampler` for use with Turing models, keeping its adaptation.
`kwargs` are passed to `externalsampler`.

```julia
sample(model, adapting_externalsampler(RobustAdaptiveMetropolis()), 10_000;
       num_warmup=5_000, check_model=false)
```
"""
adapting_externalsampler(sampler; kwargs...) =
    AdaptingExternalSampler(externalsampler(sampler; kwargs...))

function AbstractMCMC.step(
    rng::AbstractRNG, model::DynamicPPL.Model, spl::AdaptingExternalSampler; kwargs...
)
    return AbstractMCMC.step(rng, model, spl.inner; kwargs...)
end

function AbstractMCMC.step(
    rng::AbstractRNG,
    model::DynamicPPL.Model,
    spl::AdaptingExternalSampler,
    state::Turing.Inference.TuringState;
    kwargs...,
)
    return AbstractMCMC.step(rng, model, spl.inner, state; kwargs...)
end

function AbstractMCMC.step_warmup(
    rng::AbstractRNG,
    model::DynamicPPL.Model,
    spl::AdaptingExternalSampler,
    state::Turing.Inference.TuringState;
    kwargs...,
)
    f = state.ldf
    _, inner_state = AbstractMCMC.step_warmup(
        rng, AbstractMCMC.LogDensityModel(f), spl.inner.sampler, state.state; kwargs...
    )
    params = AbstractMCMC.getparams(f.model, inner_state)
    transition = DynamicPPL.ParamsWithStats(
        params, f, AbstractMCMC.getstats(inner_state)
    )
    return transition, Turing.Inference.TuringState(inner_state, params, f)
end
