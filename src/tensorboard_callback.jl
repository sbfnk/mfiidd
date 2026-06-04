#=
TensorBoard streaming callback for Turing fits
Optional backend used by reference/live_monitoring.qmd

TensorBoardLogger.jl is NOT a course dependency, so this file is opt-in:
install it once with `] add TensorBoardLogger`, then

    using TensorBoardLogger
    include(joinpath(@__DIR__, "..", "src", "tensorboard_callback.jl"))

This mirrors the optional TensorBoard extension in
epiforecasts/BVDOutbreakSize, where the heavy dependency only loads when
you ask for it. Pass the returned callback to `sample`:

    cb = tensorboard_callback("logs/run1")
    sample(model, NUTS(), 1000; callback = cb)

then view the run with `tensorboard --logdir logs/run1`.
=#

# Gate on TensorBoardLogger so an informative error is raised when it is
# not installed, rather than a bare `UndefVarError` on `log_value`.
try
    using TensorBoardLogger: TBLogger, log_value, log_histogram
catch
    error(
        "tensorboard_callback needs TensorBoardLogger.jl. Install it with " *
        "`] add TensorBoardLogger` and load it with `using TensorBoardLogger` " *
        "before including this file."
    )
end

"""
    tensorboard_callback(logdir; every = 1, histograms = true)

Build a per-iteration `callback` for `sample` that streams a Turing fit to
TensorBoard. On every kept draw it reads the sampler `transition` and logs
each real-valued quantity to the `TBLogger` at `logdir` under two grouped
tag prefixes so the dashboard stays navigable:

- `params/<name>` — every sampled parameter
- `diagnostics/<name>` — log-density (`logjoint`), divergence flag
  (`numerical_error`), step size, acceptance rate, tree depth, ...

Each scalar streams every step as a `.../value` time series. With
`histograms = true` (the default) a running histogram of the draws so far
is also logged every `every` steps as `.../distribution`, populating the
TensorBoard HISTOGRAMS and DISTRIBUTIONS dashboards. Set
`histograms = false` for scalar traces only, or raise `every` to log
histograms less often.

The whole body is wrapped in `try`/`catch` so a transition that does not
expose these fields can never abort the fit. A single logger is shared
across the threads `MCMCThreads()` spawns; writes are guarded by a
`ReentrantLock`. Use one chain for clean live traces — parallel chains
share the logger and interleave.
"""
function tensorboard_callback(
        logdir::AbstractString; every::Integer = 1, histograms::Bool = true)
    logger = TBLogger(logdir)
    lk = ReentrantLock()
    history = Dict{String, Vector{Float64}}()
    return function (rng, model, sampler, transition, state, iteration;
            kwargs...)
        try
            Base.@lock lk begin
                _emit!(logger, history, "params/", pairs(transition.params),
                    iteration, every, histograms)
                _emit!(logger, history, "diagnostics/", pairs(transition.stats),
                    iteration, every, histograms)
            end
        catch
            # A streaming callback must never bring down a fit.
        end
        return nothing
    end
end

# Log each real-valued entry of `entries` as a `<prefix><name>/value` scalar
# every step, and, when `histograms`, a `<prefix><name>/distribution`
# histogram of the accumulated draws every `every` steps.
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
