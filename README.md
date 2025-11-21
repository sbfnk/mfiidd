
# Model fitting and inference for infectious disease dynamics

This repository contains course materials for teaching Bayesian inference methods for epidemiological models at the London School of Hygiene & Tropical Medicine.

## About this course

This is a modernized version of the original [MFIIDD course](https://sbfnk.github.io/mfiidd/), reimagined for contemporary statistical practice. The course teaches both the theory and implementation of model fitting and inference using production-ready tools.

### Why the transition?

The original R-based course relied on `fitR`, a custom package developed for teaching purposes. While pedagogically useful, this approach had significant drawbacks:

- **Outdated methods**: Inference techniques have advanced substantially in recent years
- **Poor maintenance**: The teaching package became difficult to maintain and update
- **Limited performance**: Implementations were not optimized for real-world use
- **Unintended consequences**: Students adopted the teaching code for actual research, despite it not being designed for production use

### Modern approach

This revised course uses **Turing.jl**, a mature probabilistic programming language, to teach the same core concepts while providing students with:

- Modern, production-quality inference tools (NUTS, HMC, particle MCMC)
- High-performance implementations suitable for real research
- Active community support and comprehensive documentation
- Skills directly transferable to research applications

Students learn the same fundamental concepts of infectious disease modeling and Bayesian inference, but now using tools they can confidently apply in their own work.

## Course structure

Course materials are organized in:
- **`sessions/`**: Quarto documents (`.qmd`) with Julia code for hands-on sessions
- **`Rmd/`**: Legacy RMarkdown materials (being phased out)
- **`src/`**: Core Julia modeling framework with Turing.jl integration

## Local development

To render the course website locally, use Quarto:

```bash
quarto render
```

This generates HTML pages in the `_site/` directory. Changes pushed to the `main` branch are automatically deployed to the course website.
