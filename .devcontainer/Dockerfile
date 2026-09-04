FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

# Julia is installed through juliaup, the same route the course asks
# participants to use locally, so the cloud environment and a laptop end up
# running the same thing.
ENV JULIAUP_DEPOT_PATH=/usr/local/juliaup \
    JULIA_DEPOT_PATH=/home/vscode/.julia \
    PATH=/usr/local/juliaup/bin:$PATH

RUN curl -fsSL https://install.julialang.org \
      | sh -s -- --yes --default-channel 1.12 --path /usr/local/juliaup \
 && chmod -R a+rx /usr/local/juliaup \
 && julia --version

# Quarto, so a participant can render a session locally if they want to see
# their own edits. Not needed to work through the practicals.
ARG QUARTO_VERSION=1.8.27
RUN curl -fsSL -o /tmp/quarto.deb \
      "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
 && apt-get update \
 && apt-get install -y --no-install-recommends /tmp/quarto.deb \
 && rm -rf /tmp/quarto.deb /var/lib/apt/lists/*

USER vscode
