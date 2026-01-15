#!/bin/bash
#SBATCH --job-name=pigeons_pmcmc
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=pigeons_%j.out
#SBATCH --error=pigeons_%j.err

# Pigeons can be slow per sample but may give better mixing than MH

# Change to project directory
cd $SLURM_SUBMIT_DIR

echo "Starting Pigeons PMMH at $(date)"
echo "Running on node: $(hostname)"
echo "Julia version: $(julia --version)"

# Run the Pigeons script
julia --project=. scripts/run_pigeons.jl

echo "Finished at $(date)"
