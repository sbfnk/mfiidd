# Contributing

Course materials for "Model fitting and inference for infectious disease dynamics" (MFIIDD).

## Repository structure

- **`sessions/`** — Quarto documents (`.qmd`) with Julia code for each teaching session
- **`sessions/slides/`** — the revealjs lecture decks that introduce each session
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

This generates HTML in `_site/`. Changes pushed to `main` are automatically
deployed; pull requests get a preview at `previews/PR-<number>/`.

Rendering the whole site takes a while because every session executes its Julia.
To iterate on one file, render just that file:

```bash
quarto render sessions/mcmc.qmd
```

Execution results are cached in `_freeze/`, so an unchanged session is not
re-run. If a session's output looks stale after a package upgrade, delete its
directory under `_freeze/` and render again.

### Checks

`pre-commit` runs the formatter and the file hygiene hooks. Install it once:

```bash
pre-commit install
```

Julia code is formatted with JuliaFormatter following `.JuliaFormatter.toml`
(4-space indent, 92-character margin). CI runs the same checks and then renders
the site, so a render error in any session will fail the build.

Package version bounds are kept up to date by CompatHelper, which opens a pull
request when a dependency releases a new version. GitHub's dependabot does not
support Julia, which is why CompatHelper is there instead.

### Working with course materials

- Models use Turing's `@model` macro with `~` for priors and likelihoods
- Sampling: `chain = sample(model, NUTS(), n_samples)`
- Variational inference: `q, stats = vi(model, q_init, n_iter)`
- Data: `CSV.read(datadir("filename.csv"), DataFrame)`
- Diagnostics: StatsPlots.jl for trace plots, density plots, convergence checks

### Time indexing and the initial state

Every session that fits the Tristan da Cunha data shares these two conventions.
A session that departs from either will disagree with the rest of the course by a day or by a few cases, and the disagreement is quiet enough to survive a review.

`times` starts at `0.0`, the day before the first observation, so the solver has a day to integrate over before the first count.
The incidence on that leading day is then discarded: write `traj.Inc[i + 1]` in a likelihood, and return `Inc[2:end]` or a bare `diff(...)` from a simulator wrapper.
Observation `i` then lines up with the day it counts.
Plot predictions against `flu_tdc.time`, because the trimmed vector is one element shorter than `times`.

The initial state is `S = 279`, `I = 2`, `R = 3`, giving N = 284, the island population.
Of the five islanders who landed on day 0, three had been ill during the eight day voyage and are past their infectious period, so they start recovered; the two who fell ill on landing start infectious.
The 312 reported cases exceed N because islanders were infected more than once, which is the observation the SEITL session is built on and the structural failure the single-wave sessions diagnose.

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
| Open-ended answer | `callout-tip title="What to look for" collapse="true"` | The answer to an exercise that has no single right answer: what a good reader should notice, and why it matters. |
| Gotcha / pitfall | `callout-warning` | "This will bite you": failure modes, computational constraints, common mistakes, stochastic noise to expect. |
| Key takeaway | `callout-important` | Concepts the student must hold for the rest of the session. Budget: ≤2 per session. |

Other rules:

- End each session with a `callout-tip title="Learning points"` box of bullet points. Do not also write a prose `# Summary` section — pick one.
- Never wrap section headings (`# Setup`, `# Introduction`, etc.) inside a callout. Callouts are visual; sections are structural. If code blocks are bulky, box the code block with its own title and leave the section heading in the document outline.
- "Coming from R?" / "Coming from Stan?" is always a collapsed blue note.

## Slide decks

Each taught session has a revealjs deck in `sessions/slides/` that the instructor
presents before the practical. Shared settings live in
`sessions/slides/_metadata.yml` and apply to every deck in the directory, so an
individual deck's front matter only needs a `title`, `subtitle` and `footer`.

Styling is in `sessions/slides/mfiidd.scss`. A few classes are worth knowing:

- `[text]{.alert}` for the one emphasised term in a sentence, as `\alert{}` did
  in the old Beamer decks
- `{.step}` on a figure that is one frame of a sequence, so successive frames
  are drawn at a consistent size
- `{.tall}` for portrait figures, `{.mini}` for a small figure above code
- `[source]{.cite}` for a small grey attribution under a figure

### Stepped figures

A sequence of figures shown in the same place goes in an `.r-stack`, with the
first frame fading out as the second appears:

```markdown
::: {.r-stack}
![](images/frame_1.svg){.step .fragment .fade-out fragment-index=1}
![](images/frame_2.svg){.step .fragment fragment-index=1}
:::
```

The first frame needs `.fade-out` rather than plain `.fragment`, or the slide
starts blank and the audience sees nothing until you press space.

The `.r-stack` rules in the SCSS include a `display: contents` on the wrapping
paragraph. Pandoc collects a run of images into a single `<p>`, which stops
reveal's own overlay rule from reaching the images; that declaration puts them
back. The file has a comment explaining it. If you change those rules, check the
result by putting two images in an `.r-stack` with no fragments — they must sit
on top of each other.

### Figures

Prefer SVG. The decks are projected, and the source figures are line drawings
that stay sharp at any size. `sessions/slides/images/` holds the originals.

### Viewing a deck

```bash
quarto preview sessions/slides/mcmc.qmd
```

Right arrow moves between top-level (`#`) sections; the space bar walks every
slide including the `##` slides nested under them. Press `s` for speaker notes
and `o` for the overview grid.

## Adding a session

A new session touches more than one file. In order:

1. Write `sessions/<name>.qmd`. Set `order:` in the front matter to place it in
   the sidebar, and open with the estimated-time banner the other sessions use.
2. Add a deck at `sessions/slides/<name>.qmd` if the session is taught, ending
   with a "Your Turn" slide and a link back to the session.
3. Add the session to the `sessions` sidebar in `_quarto.yml`.
4. Add it to the course map in `reference/structure.qmd`, including its
   prerequisites.
5. If it is taught rather than self-study, give it a slot in
   `reference/sessions.qmd` with per-part timings that add up to the slot.
6. Render locally before opening the pull request. The session must execute from
   a clean `_freeze/`.

## Exercises

Sessions carry their own exercises. Look at `sessions/model_checking.qmd` for the
pattern: a `callout-tip` posing the task, followed by a collapsed answer. Aim for
two or more per session.

**Exercises are conceptual, and all the Julia is given complete.** There are
deliberately no fill-in-the-blank coding tasks, unlike the older R version of the
course. Students were spending practical time debugging code rather than learning
about inference, and the debugging taught them nothing about the subject. Write
exercises that cannot break at the coding stage: predict the output before running
the cell, interpret a diagnostic, diagnose why a fit failed, choose between
options and justify it.

**Every exercise needs a collapsed answer.** Most of the course is worked through
by people reading on their own, and an exercise with no answer gives them a
question and silence. Use `title="Solution"` where there is a determinate answer
and `title="What to look for"` where it is a matter of judgement — calling a
judgement call a solution overstates it.

Answers should be specific to the session's own numbers. "The posterior is wider"
is worth little; "the stochastic model gives wider posteriors because it has a
second source of variation available, so a narrower posterior is not a better one
when the narrowness comes from a mechanism you left out" is worth reading.

## Data and long-running computations

Anything that takes more than a couple of minutes is precomputed. The scripts
that generate those results live in `scripts/`, and the sessions load the saved
output from `data/`. If you change a model that a precomputed chain depends on,
re-run the script and commit the new `.rdata` file alongside the change.
