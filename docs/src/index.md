# MFIIDD

API reference for the `MFIIDD` Julia package, the supporting library for the
[Model fitting and inference for infectious disease dynamics](https://sbfnk.github.io/mfiidd/)
course.

The course session notebooks under `sessions/` are rendered by Quarto and live
at the root of the course site. This Documenter site is published as a
sub-path (`/api`) and is the authoritative reference for the package's
exported functions and types.

See the [API](api.md) page for the full reference.

## Building locally

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
```

The output is written to `docs/build/`.
