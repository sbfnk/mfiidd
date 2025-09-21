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
The course uses a modern Bayesian approach for epidemiological models:
- **`FitModel` struct**: Contains name, state variables, parameters, `Distribution` objects for priors, simulation and likelihood functions
- **Turing.jl integration**: Direct translation to `@model` functions for MCMC sampling
- **Distribution objects**: Elegant prior specification using `Distributions.jl` (e.g., `Uniform(1.0, 100.0)`, `Normal(2.5, 0.5)`)
- **Likelihood patterns**: Using `~` syntax compatible with Turing.jl (e.g., `obs[i] ~ Poisson(expected_cases[i])`)
- **MCMC sampling**: Built-in support for NUTS, HMC, and other modern samplers

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
- **src/FitModels.jl**: Core epidemiological modeling framework with Turing.jl integration
  - `FitModel` struct with `Distribution` objects for priors
  - `@model` functions for Turing.jl compatibility
  - Functions: `log_prior()`, `log_likelihood()`, `log_posterior()`, `sample_posterior()`
  - Visualization: `plot_trajectory()`, `plot_fit()`, `plot_posterior_predictive()`
- **src/SampleData.jl**: Sample epidemic datasets and model loading functions
- **qmd/**: Quarto documents demonstrating Bayesian workflows with Julia
- **test_julia_code.jl**: Comprehensive tests including Turing.jl integration

### Working with Course Materials
- **Modern Bayesian workflow**: Prior specification → MCMC sampling → Posterior analysis
- **Elegant syntax**: Use `Distribution` objects directly (e.g., `R_0 ~ Uniform(1.0, 100.0)`)
- **MCMC sampling**: `chain = sample_posterior(model, data, init_state, n_samples=1000)`
- **Diagnostics**: Use `StatsPlots.jl` for trace plots, density plots, and convergence checks
- **Posterior predictive checks**: Built-in functions for model validation
- Course data loaded via `load_models()` and `load_epi_data()` functions
- Use `include("../src/FitModels.jl")` to load the full Bayesian modeling framework

## Development Notes

- The repository is in active transition from R to Julia-based instruction
- Preserve existing R code structure when adding Julia equivalents
- Course materials assume familiarity with epidemiological modeling concepts
- Models follow a standardized structure for educational consistency
- Cache directories and generated HTML should not be committed