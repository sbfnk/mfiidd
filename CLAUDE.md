# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains course materials for "Model fitting and inference for infectious disease dynamics" (MFIIDD) at the London School of Hygiene & Tropical Medicine. The project is transitioning from R to Julia-based instruction while maintaining R-based examples and legacy code.

## Architecture

### Content Structure
- **Rmd/**: RMarkdown source files for course content and exercises
  - Main course pages: `index.Rmd`, `introduction.Rmd`, `mcmc.Rmd`, `mcmc_diagnostics.Rmd`, etc.
  - `_site.yml`: Website navigation and configuration
  - Course exercises and solutions with pattern `*_example.Rmd` and `*_solution.Rmd`
  - Projects in `Projects/` subdirectory
- **scripts/**: R code organized by purpose
  - `models/`: Epidemiological model definitions (SIR, SEITL variants)
  - `generate-data/`: Data generation scripts for exercises
  - `snippets/`: Reusable code chunks for course materials
  - `figures/`: Figure generation scripts
- **docs/**: Documentation including Julia reference (`rintro.org`)
- **rds/**: Serialized R data objects for course exercises
- **MCMC slides/**: LaTeX beamer slides and assets

### Technology Stack
- **Primary Language**: Julia (v1.11.7+) with modern Bayesian inference
- **Legacy Content**: R with RMarkdown
- **Core Dependencies**: 
  - **Turing.jl**: Probabilistic programming and MCMC sampling
  - **Distributions.jl**: Elegant prior specification and likelihood calculation
  - **DifferentialEquations.jl**: ODE solving for epidemiological models
  - **Plots.jl/StatsPlots.jl**: Visualization and MCMC diagnostics
  - **MCMCChains.jl**: MCMC chain analysis and diagnostics
- **Website Generation**: Quarto for Julia code blocks, RMarkdown for legacy content

### Model Framework
The course uses Turing.jl for Bayesian inference:
- **Turing.jl `@model` functions**: Define models with priors and likelihoods inline
- **Distribution objects**: Prior specification using `Distributions.jl` (e.g., `Uniform(1.0, 100.0)`, `Normal(2.5, 0.5)`)
- **Likelihood patterns**: Using `~` syntax (e.g., `obs[i] ~ Poisson(expected_cases[i])`)
- **MCMC sampling**: Built-in support for NUTS, HMC, and other modern samplers
- **Variational inference**: Fast approximate inference with `Turing.Variational`

## Common Development Commands

### Build and Test
```r
# Generate website locally (legacy R content)
rmarkdown::render_site("Rmd/")
```

```julia
# Test Julia code translation
julia --project=. test_julia_code.jl

# Activate project environment
julia --project=.

# Install/update dependencies
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

### Julia Setup
- Ensure Julia version ≥ 1.11.7
- Use `Project.toml` and `Manifest.toml` for dependency management
- Key dependencies: DifferentialEquations, Distributions, DataFrames, Plots
- Development dependencies include LanguageServer and Revise for interactive development

### Julia Code Structure
- **sessions/**: Quarto documents (`.qmd`) with Julia code blocks for each teaching session
- **data/**: CSV data files (e.g., `flu_tdc_1971.csv` for Tristan da Cunha outbreak)
- **_freeze/**: Cached execution results for Quarto rendering

### Working with Course Materials
- **Bayesian workflow**: Prior specification → MCMC sampling → Posterior analysis
- **Syntax**: Use Turing's `@model` macro with `~` for priors and likelihoods
- **MCMC sampling**: `chain = sample(model, NUTS(), n_samples)`
- **VI sampling**: `q, stats = vi(model, q_init, n_iter)`
- **Diagnostics**: Use `StatsPlots.jl` for trace plots, density plots, and convergence checks
- **Data**: Load with `CSV.read(datadir("filename.csv"), DataFrame)`

## Development Notes

- The repository is in active transition from R to Julia-based instruction
- Preserve existing R code structure when adding Julia equivalents
- Course materials assume familiarity with epidemiological modeling concepts
- Models follow a standardized structure for educational consistency
- Cache directories and generated HTML should not be committed