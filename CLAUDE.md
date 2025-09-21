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
- **Primary Language**: Transitioning to Julia (v1.11.7+)
- **Legacy Content**: R with RMarkdown
- **Dependencies**: Julia packages listed in `Project.toml` (DifferentialEquations, LanguageServer, Revise, SymbolServer)
- **Website Generation**: RMarkdown site using `rmarkdown::render_site()`

### Model Framework
The course uses a structured approach for epidemiological models with these components:
- `fitmodel` objects containing: name, state variables, parameters, simulation functions, prior distributions, likelihood functions
- Deterministic and stochastic variants of each model
- Standardized interface for MCMC, ABC, and particle MCMC methods

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
- **src/FitModels.jl**: Core epidemiological modeling framework with `FitModel` struct
- **src/SampleData.jl**: Sample epidemic datasets and model loading functions
- **qmd/**: Quarto documents with Julia code blocks (converted from Rmd)
- Models follow pattern: simulation functions, prior/likelihood evaluation, observation generation

### Working with Course Materials
- RMarkdown files use `cache=TRUE` for expensive computations (legacy)
- Quarto files use Julia execution environment
- Course data loaded via `load_models()` and `load_epi_data()` functions
- Exercise solutions follow naming pattern: `*_example.qmd` → `*_solution.qmd`
- Use `include("../src/FitModels.jl")` to load Julia modeling framework

## Development Notes

- The repository is in active transition from R to Julia-based instruction
- Preserve existing R code structure when adding Julia equivalents
- Course materials assume familiarity with epidemiological modeling concepts
- Models follow a standardized structure for educational consistency
- Cache directories and generated HTML should not be committed