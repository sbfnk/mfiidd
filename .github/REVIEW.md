# Reviewing a change to mfiidd

What to look for in a change to the MFIIDD course materials specifically — a
Quarto site of teaching sessions (`sessions/*.qmd`) whose Julia executes at
render time, using Turing.jl for the inference the course teaches. The reviewing
method — how to scope a review, what counts as a finding, how to report and
suggest the fix, the trust rules — is the org half of this spec
(`mfiidd/.github` → `REVIEW.md`), and a review follows both. The repository's
conventions are in `CONTRIBUTING.md` and `CLAUDE.md`.

**This is teaching material, and that changes the bar.** Code here is written to
be read by someone learning the method, not to be idiomatic or fast. A loop
written out longhand where a comprehension would do, a model written inline
rather than factored into a function, a variable named for the notation in the
lecture rather than for Julia style — these are pedagogy, not defects. Do not
suggest a refactor whose only benefit is elegance. The question is whether a
participant following the session would be misled, not whether an experienced
Julia programmer would have written it differently.

## What to look for

- **The statistics being taught**: a prior that does not match what the text says
  it is, a likelihood that does not correspond to the stated observation model, a
  `@model` whose parameters are sampled but never used. An error here teaches the
  wrong thing, which is worse than a bug that merely breaks.
- **Session narrative and code out of step**: prose describing a parameter,
  result, figure or number that the code no longer produces. These drift apart
  whenever a model is retuned, and the text is what participants trust.
- **Reproducibility of sampling**: MCMC, particle filters and ABC all draw
  randomly. A chain, figure or quoted result that changes between renders needs a
  seed; a quoted diagnostic (`ess`, `rhat`, an acceptance rate) that was true of
  one run and is now hard-coded in the prose will go stale silently.
- **Render cost**: sessions execute on every render, so a change that adds
  iterations, particles or a wider grid makes the site slower for everyone.
  Long-running computation belongs in `scripts/` with its output committed, not
  inline in a session.
- **Exercises left deliberately incomplete**: a gap a participant is meant to fill
  is not a bug. Check the surrounding text before flagging missing code.
- **Cross-session consistency**: the SEITL model, its parameter names and the
  datasets in `data/` recur across sessions. A change to one session's definition
  that leaves the others alone breaks the thread the course follows.
- **Things that must not be committed**: `_freeze/`, `_site/`, generated HTML.
