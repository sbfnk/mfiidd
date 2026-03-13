
# Model fitting and inference for infectious disease dynamics

Course materials for teaching Bayesian inference methods for infectious disease models, developed at the London School of Hygiene & Tropical Medicine. The course uses Julia with the [Turing.jl](https://turinglang.org/) probabilistic programming framework.

## What the course covers

The course bridges the gap between statistical inference methods and infectious disease modelling. It covers:

- Fitting deterministic and stochastic compartmental models to data
- MCMC sampling (NUTS, HMC) and convergence diagnostics
- Model checking with prior and posterior predictive checks
- Particle filters and particle MCMC for stochastic models
- Observation models and Approximate Bayesian Computation
- Variational inference and universal differential equations (additional sessions)

No prior Julia experience is required. The materials include guidance for those coming from R throughout.

## Repository structure

- **`sessions/`**: Quarto documents (`.qmd`) with Julia code for each teaching session
- **`data/`**: Datasets used in practicals
- **`Rmd/`**: Legacy RMarkdown materials from the [previous version of the course](https://sbfnk.github.io/mfiidd.archive/)

## Local development

To render the course website locally, use [Quarto](https://quarto.org/):

```bash
quarto render
```

This generates HTML pages in the `_site/` directory. Changes pushed to the `main` branch are automatically deployed to the course website.

## Licence

[MIT](LICENSE)
