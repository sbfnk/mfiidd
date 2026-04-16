# Contributing

Course materials for "Model fitting and inference for infectious disease dynamics" (MFIIDD) at the London School of Hygiene & Tropical Medicine.

## Repository structure

- **`sessions/`** — Quarto documents (`.qmd`) with Julia code for each teaching session
- **`data/`** — CSV datasets used in practicals
- **`figures/`** — Diagrams used in sessions (compartmental model figures, etc.)
- **`scripts/`** — Julia scripts for long-running computations (e.g., PMCMC chain generation)
- **`src/`** — Shared Julia source code
- **`reference/`** — Reference materials
- **`_freeze/`** — Cached Quarto execution results (not committed)

## Technology stack

- **Julia** (v1.11.7+) with Turing.jl for probabilistic programming
- **Quarto** for rendering the course website
- Key dependencies: DifferentialEquations.jl, Distributions.jl, Plots.jl/StatsPlots.jl, MCMCChains.jl

## Development

### Setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Render the site locally

```bash
quarto render
```

This generates HTML in `_site/`. Changes pushed to `main` are automatically deployed.

### Working with course materials

- Models use Turing's `@model` macro with `~` for priors and likelihoods
- Sampling: `chain = sample(model, NUTS(), n_samples)`
- Variational inference: `q, stats = vi(model, q_init, n_iter)`
- Data: `CSV.read(datadir("filename.csv"), DataFrame)`
- Diagnostics: StatsPlots.jl for trace plots, density plots, convergence checks

## Style

- Use British English in all written content
- Models follow a standardised structure for educational consistency
- Cache directories and generated HTML should not be committed
