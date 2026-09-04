# Codespaces environment

This gives a participant a browser-based VS Code with Julia and every course
package already installed, for anyone whose local installation will not
cooperate, or who would rather not download a gigabyte in a teaching room.

## Why the commands sit where they do

Instantiating the project downloads around 1 GB and precompiles 774
dependencies, measured at 28 minutes on the two-core machine declared below. That work happens in `onCreateCommand`,
which a **prebuild** runs and captures in the image, so opening a codespace
takes seconds instead. Moving it to `postCreateCommand` would lose the saving
entirely, because that runs when a participant opens their codespace rather
than when the image is built.

`updateContentCommand` just loads `MFIIDD` once. It is cheap and it fails the
prebuild loudly if the image is broken, rather than leaving someone to discover
it on the morning.

A prebuild is configured on `main`, so that work is already done and a codespace
opens without it. Prebuilds require the owning organisation to be on GitHub Team,
which `mfiidd` reached through a GitHub Education faculty coupon, and they need
the organisation's Codespaces spending limit set above zero, because the prebuilt
image is stored and storage is billed.

Enable or change it under **Settings → Codespaces → Set up prebuild**, on `main`,
triggered by configuration change. The first build took 72 minutes.

Two cases still pay the full 28 minutes, so it is worth knowing they exist: a
codespace opened on a branch rather than `main`, and one opened after
`Project.toml`, `Manifest.toml` or this directory changed but before the rebuild
finishes. Neither fails; they are just slow.

## What opens first

`customizations.codespaces.openFiles` opens the course overview and the first
session instead of leaving the participant on `README.md`, which is written for
people working on the course rather than taking it. It applies when the codespace
is created, not on every reconnect.

## Machine size

Two cores and 8 GB, which is `basicLinux32gb`. With the prebuild in place there is
no setup wait to shorten, so the only thing more cores would change is the rate at
which a participant spends their allowance.

A codespace on a public repository is billed against the user's **personal** free
allowance of 120 core-hours a month. Four days of teaching at six and a half hours
a day costs about 52 core-hours on two cores and about 104 on four.

Note that none of this is billed to the organisation. Org billing for
codespaces requires Team or Enterprise Cloud, so on the free plan every
codespace is charged to the personal allowance of whoever created it. There is
also no way to create one on someone else's behalf: a codespace belongs to the
account that made it.

## Comments

`devcontainer.json` permits comments, but this repository's `check-json`
pre-commit hook parses it as strict JSON and rejects them, which is why the
reasoning lives here instead of inline.

## Versions

The Dockerfile pins the juliaup channel to 1.12, matching the `julia = "1.12"`
compat bound in `Project.toml`. If that bound moves, change it here too.

Quarto is deliberately **not** pinned. `build-and-deploy.yml` uses
`quarto-actions/setup` with no version, so the published site is built with
whatever is current, and pinning here would only guarantee that a session
rendered in a codespace differs from the site it is being compared against.
Both now track the latest release.
