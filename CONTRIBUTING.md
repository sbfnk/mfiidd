# Contributing

Course materials for "Model fitting and inference for infectious disease dynamics" (MFIIDD).

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

## Callout conventions

Callouts are signals, not decoration. Use prose by default; reach for a callout only if the content would otherwise be missed or skipped. One colour per purpose:

| Role | Syntax | Use for |
|------|--------|---------|
| Metadata banner | `callout-note appearance="simple"` | Session-level metadata only (e.g. estimated time). One per session, at the top. |
| Optional context | `callout-note collapse="true"` | Background a linear reader can skip: "Coming from R/Stan", derivations, code deep-dives, linked implementation files. |
| Exercise / Hint / Solution | `callout-tip` (`collapse="true"` for hints and solutions) | Anything the student *does* rather than reads. |
| Gotcha / pitfall | `callout-warning` | "This will bite you": failure modes, computational constraints, common mistakes, stochastic noise to expect. |
| Key takeaway | `callout-important` | Concepts the student must hold for the rest of the session. Budget: ≤2 per session. |

Other rules:

- End each session with a `callout-tip title="Learning points"` box of bullet points. Do not also write a prose `# Summary` section — pick one.
- Never wrap section headings (`# Setup`, `# Introduction`, etc.) inside a callout. Callouts are visual; sections are structural. If code blocks are bulky, box the code block with its own title and leave the section heading in the document outline.
- "Coming from R?" / "Coming from Stan?" is always a collapsed blue note.
