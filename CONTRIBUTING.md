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

Use these two conventions in any session that fits the Tristan da Cunha data.
A session that departs from either disagrees with the others by a day or by a few cases, quietly enough to survive a review, so it is worth checking a session against this section rather than against whichever session happens to be open.

`times` starts at `0.0`, one point before the first observation, because the incidence over observation day 1 is `S(0) - S(1)` and there is no way to compute it without a day 0.
That leading day is then dropped, in exactly one of two places.
Either the simulator wrapper trims it, returning `Inc[2:end]` or a bare `diff(...)`, and the likelihood indexes `inc[i]`; or the wrapper returns the full vector and the likelihood indexes `traj.Inc[i + 1]`.
Do one or the other, never both: a wrapper that trims and a likelihood that also shifts ask for one index past the end of the vector, and the session fails to render with a `BoundsError`.
Either way observation `i` lines up with the day it counts.
Plot predictions against `flu_tdc.time`, or against `collect(times)[2:end]` in a cell that runs before the data are loaded, because the trimmed vector is one element shorter than `times`.

The initial state is `S = 279`, `I = 2`, `R = 3`, giving N = 284, the island population.
It is the island as the ship lands, at the start of day 1: `flu_tdc_1971.csv` dates `time = 1` as 13 August, the day the ship arrived.
Of the five islanders it brought back, three had been ill during the eight day voyage and are past their infectious period, so they start recovered; the two who fell ill on landing start infectious.
The 312 reported cases exceed N because islanders were infected more than once, which is the observation the SEITL session is built on and the structural failure the single-wave sessions diagnose.

## House style

### Prose

- **British English throughout**: modelling, behaviour, summarise, visualise,
  centre. This applies to prose, comments, commit messages and the text of figure
  labels. It does not apply to code: identifiers, CSS properties and Quarto or
  reveal attribute values keep their own spelling, so `fig-align: center`,
  `text-align: center` and the `{.center}` class stay as they are.
- **No contractions.** Write "do not", "we will", "it is", "let us". The course
  reads as written teaching material rather than as a transcript of the lecture.
  Only `julia_ecosystem.qmd` currently follows this throughout; the other eleven
  sessions carry between 4 and 37 contractions each, so applying the rule to the
  existing material is a substantial sweep rather than a tidy-up.
- **One sentence per line** in the `.qmd` source. Start each sentence on a new
  line and do not wrap within a sentence, however long it runs. Markdown joins
  the lines when rendering, so the output is unchanged, and a reworded sentence
  then shows up as a single changed line instead of a reflowed paragraph. The
  gain is entirely in review: prose diffs become readable.
  - A sentence beginning with a numeral and a full stop ("1998. That year...")
    is parsed as an ordered-list item, but only when it is the first line of a
    paragraph. Mid-paragraph it stays inline, which is where most sentences end
    up under this rule. Reword only when it bites.
  - Tables, YAML, callout delimiters and code blocks are exempt.
- **Sentence case for headings**: "Prior predictive checks", not "Prior
  Predictive Checks". Proper nouns keep their capitals ("Sampling with
  DynamicHMC.jl").
- **"We" for what the session does, "you" for what the reader does.** "We fit
  the model with NUTS" describes the worked example; "you should see a warning"
  addresses the reader. Both are used throughout and neither should be purged.
- Avoid "X, not Y" and "X, rather than Y" unless a reader who has not seen the
  drafting history would recognise Y as a real alternative worth weighing. State
  the positive claim on its own.
- **End-of-line comments in Julia use `##`**, as in `using Random ## for random
  numbers`. This is a course convention rather than a Julia one, which uses a
  single `#`, but the sessions are near-unanimous: `##` outnumbers `#` by about
  sixty to four. Follow it, or change it everywhere in one sweep and update this
  line. Do not fix it session by session, which only splits the corpus further.

### Session structure

Every session has the first six of these, in this order. The closing sections
vary across the corpus, and the ordering given here is the target rather than a
description:

1. YAML front matter. Every session carries `title`, `format: html`,
   `engine: julia` and `order`; without the last two the Julia does not execute
   and the rendered page is codeless. `subtitle` is optional and four sessions
   use it.
2. The estimated-time banner, a `callout-note appearance="simple"`. Every session
   gives an estimated time. Add a level marker where one applies, as
   `**Advanced session**` or `**Optional session**`, and `**Requires**:` where
   the session genuinely depends on earlier ones, which is five of the twelve
   today. Linking the prerequisites is welcome but not the norm yet.
3. `# Introduction`, motivating the session from where the previous one stopped.
4. `# Objectives`, a numbered list of what the reader will be able to do.
5. `# Setup`, loading packages and data. Data loading belongs here rather than
   in a top-level section of its own.
6. The body.
7. A `callout-tip title="Learning points"` box of bullet points, once per
   session, at the end of the body. The Callout conventions section below says the
   same thing; this is the only place the placement is specified.
8. Optionally `# Going further` and `# Next session`.
9. `# References` last, holding the `::: {#refs}` block.

Four sessions do not match points 8 and 9 today: `observation_models.qmd` has no
`# References` section at all, `abc.qmd` ends `Going further` → `References` →
`Next session`, `seitl.qmd` ends `References` → `Going further`, and
`universal_differential_equations.qmd` places `Going further` before its Learning
points box. Bring a session into line when you are editing it for another reason
rather than making a sweep of its own.

### Exercise placement

- **Interleave by default**, placing each exercise immediately after the material
  it tests. Most sessions do this and it keeps the exercise next to the code it
  refers to.
- **Advanced sessions may collect exercises** in a single `# Exercises` section
  before `# Going further`, which is what the variational inference and universal
  differential equations sessions do. Pick one arrangement per session rather
  than splitting between the two.
- Two or more per session, each with a collapsed answer. See "Exercises" below
  for what makes a good one.

### "Coming from R?"

Write one where R genuinely does the thing differently and the difference would
otherwise trip up an R-speaking reader: a different default, a missing library, a
convention that does not transfer. Do not add one to a session for the sake of
having one, and do not write one that only says the syntax differs. Sessions
without a real contrast are better without the note.

### Other conventions

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

## Reviewing a change

`.github/REVIEW.md` says what a review of this repository looks for, and what is
not worth reporting. Follow it yourself, or point a coding assistant at it.
