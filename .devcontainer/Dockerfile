FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

# Julia is installed through juliaup, the same route the course asks
# participants to use locally, so the cloud environment and a laptop end up
# running the same thing.
ENV JULIAUP_DEPOT_PATH=/usr/local/juliaup \
    JULIA_DEPOT_PATH=/home/vscode/.julia \
    PATH=/usr/local/juliaup/bin:$PATH

RUN curl -fsSL https://install.julialang.org \
      | sh -s -- --yes --default-channel 1.12 --path /usr/local/juliaup \
 && chown -R vscode:vscode /usr/local/juliaup \
 && julia --version

USER vscode
