#!/bin/bash
#SBATCH --job-name=pmcmc_chains
#SBATCH --time=48:00:00   # ~15 hours per model, run back to back
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=pmcmc_%j.out
#SBATCH --error=pmcmc_%j.err

# Julia is installed via juliaup on the cluster, so it is already on PATH

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
