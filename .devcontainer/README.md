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

There is no prebuild, and there cannot be one. Prebuilds require the owning
organisation to be on GitHub Team or Enterprise Cloud, and `mfiidd` is on the
free plan, which is why the repository has no **Settings → Codespaces** page at
all.

So creating a codespace blocks on `onCreateCommand`, showing a window that looks
like it has hung. Assume this happens during the setup slot on the morning
rather than in advance: people arrive and start then, whatever the instructions
say. That is why the machine size below is chosen around the setup time and not
around the running cost.

## Machine size

Two cores and 8 GB, which is `basicLinux32gb`. Four cores would roughly halve the
28-minute setup, and it is tempting for exactly that reason, but a codespace on a
public repository is billed against the user's **personal** free allowance of 120
core-hours a month. Four days of teaching at six and a half hours a day costs
about 52 core-hours on two cores and about 104 on four, which leaves nothing for
setting up beforehand or for the rest of the month.

Creating the codespace the day before buys back the same 28 minutes and costs
nothing, so the cores are not worth spending.

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
