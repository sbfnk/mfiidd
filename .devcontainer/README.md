# Codespaces environment

This gives a participant a browser-based VS Code with Julia and every course
package already installed, for anyone whose local installation will not
cooperate, or who would rather not download a gigabyte in a teaching room.

## Why the commands sit where they do

Instantiating the project downloads around 1 GB and precompiles 491 packages,
which takes the best part of an hour. That work happens in `onCreateCommand`,
which a **prebuild** runs and captures in the image, so opening a codespace
takes seconds instead. Moving it to `postCreateCommand` would lose the saving
entirely, because that runs when a participant opens their codespace rather
than when the image is built.

`updateContentCommand` just loads `MFIIDD` once. It is cheap and it fails the
prebuild loudly if the image is broken, rather than leaving someone to discover
it on the morning.

Prebuilds are not configured in this file, and cannot be: they are a repository
setting. Enable them per branch under **Settings → Codespaces → Set up prebuild**,
and rebuild after any change to `Project.toml`, `Manifest.toml` or this directory.

**Do that before pointing anyone at a codespace.** Without a prebuild, creating one
blocks on `onCreateCommand` while it downloads the packages and precompiles them,
which was measured at 26 minutes on the 2-core machine declared here, on top of
around 2.5 GB of downloads. The participant sees a window that appears to hang,
with no indication why, which is worse than the local installation this option
exists to rescue.

## Machine size

Two cores and 8 GB is enough for the practicals. The usual reason to want more
cores is precompilation, and the prebuild has already done that.

## Comments

`devcontainer.json` permits comments, but this repository's `check-json`
pre-commit hook parses it as strict JSON and rejects them, which is why the
reasoning lives here instead of inline.

## Pinned versions

The Dockerfile pins the juliaup channel to 1.12, matching the `julia = "1.12"`
compat bound in `Project.toml`. If that bound moves, change it here too.

It also pins Quarto through the `QUARTO_VERSION` build argument, which
`build-and-deploy.yml` does not: `quarto-actions/setup` takes no version and so
tracks the latest release. The two therefore drift apart, and a session rendered
in a codespace can differ from the published site. Bump the argument when the
gap starts to matter.
