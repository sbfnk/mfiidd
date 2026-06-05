#=
TensorBoard streaming callback for Turing fits
Used by reference/live_monitoring.qmd

Logs each sampled parameter and each sampler diagnostic to TensorBoard, so a
fit (or a finished run) can be explored in the browser dashboard. Mirrors the
optional TensorBoard backend in epiforecasts/BVDOutbreakSize.

Usage:

    using TensorBoardLogger
    include(joinpath(@__DIR__, "..", "src", "tensorboard_callback.jl"))

    cb = tensorboard_callback("logs/run1")
    sample(model, NUTS(), 1000; callback = cb)

then view the run with `tensorboard --logdir logs/run1`.
=#

using AbstractMCMC: ParamsWithStats
using TensorBoardLogger: TBLogger, log_value, log_histogram

"""
    tensorboard_callback(logdir; every = 1, histograms = true)

Build a per-iteration `callback` for `sample` that streams a Turing fit to
TensorBoard. Parameters and step statistics are read through the
sampler-agnostic `AbstractMCMC.ParamsWithStats(model, sampler, transition,
state; params, stats)` interface (so the callback tracks the transition format
rather than a fixed field layout) and logged under two grouped tag prefixes so
the dashboard stays navigable:

- `params/<name>` — every sampled parameter
- `diagnostics/<name>` — log-density (`logjoint`), divergence flag
  (`numerical_error`), step size, acceptance rate, tree depth, ...

Each scalar streams every step as a `.../value` time series (the SCALARS
dashboard). With `histograms = true` (the default) a running histogram of the
draws so far is also logged every `every` steps as `.../distribution`,
populating the HISTOGRAMS and DISTRIBUTIONS dashboards. Set `histograms = false`
for scalar traces only, or raise `every` to log histograms less often.

The whole body is wrapped in `try`/`catch` so a transition that does not expose
these fields can never abort the fit. A single logger is shared across the
threads `MCMCThreads()` spawns; writes are guarded by a `ReentrantLock`. To
compare several chains, give each its own `logdir` and point TensorBoard at the
parent directory.
"""
function tensorboard_callback(
        logdir::AbstractString; every::Integer = 1, histograms::Bool = true)
    logger = TBLogger(logdir)
    lk = ReentrantLock()
    history = Dict{String, Vector{Float64}}()
    return function (rng, model, sampler, transition, state, iteration;
            kwargs...)
        try
            pws = ParamsWithStats(
                model, sampler, transition, state; params = true, stats = true)
            Base.@lock lk begin
                _emit!(logger, history, "params/", pairs(pws.params),
                    iteration, every, histograms)
                _emit!(logger, history, "diagnostics/", pairs(pws.stats),
                    iteration, every, histograms)
            end
        catch
            # A streaming callback must never bring down a fit.
        end
        return nothing
    end
end

# Log each real-valued entry of `entries` as a `<prefix><name>/value` scalar
# every step, and, when `histograms`, a `<prefix><name>/distribution` histogram
# of the accumulated draws every `every` steps.
function _emit!(logger, history, prefix, entries, iteration, every, histograms)
    for (name, value) in entries
        (value isa Real && isfinite(value)) || continue
        v = Float64(value)
        tag = string(prefix, name)
        log_value(logger, string(tag, "/value"), v; step = iteration)
        if histograms
            draws = get!(() -> Float64[], history, tag)
            push!(draws, v)
            if iteration % every == 0 && length(draws) > 1
                log_histogram(logger, string(tag, "/distribution"), draws;
                    step = iteration)
            end
        end
    end
    return nothing
end
