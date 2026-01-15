#!/bin/bash
#SBATCH --job-name=pmcmc_chains
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=pmcmc_%j.out
#SBATCH --error=pmcmc_%j.err

# Load Julia (adjust path if using module or different install location)
export PATH="$HOME/julia-1.11.2/bin:$PATH"
# Or if using juliaup:
# export PATH="$HOME/.juliaup/bin:$PATH"

# Change to project directory
cd $SLURM_SUBMIT_DIR
# Or use absolute path:
# cd /path/to/mfiidd.julia

echo "Starting PMMH chain generation at $(date)"
echo "Running on node: $(hostname)"
echo "Julia version: $(julia --version)"

# Run the chain generation script
julia --project=. scripts/generate_pmcmc_chains.jl

echo "Finished at $(date)"
