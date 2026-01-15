#!/bin/bash
#SBATCH --job-name=adaptive_pmmh
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=adaptive_pmmh_%j.out
#SBATCH --error=adaptive_pmmh_%j.err

# Adaptive MCMC learns the posterior covariance structure
# Should give much better mixing than independent proposals

cd $SLURM_SUBMIT_DIR

echo "Starting Adaptive PMMH at $(date)"
echo "Running on node: $(hostname)"
echo "Julia version: $(julia --version)"

julia --project=. scripts/adaptive_pmmh.jl

echo "Finished at $(date)"
