# Generate the pre-computed SEITL PMMH chain loaded by sessions/pmcmc.qmd.
# Run with: julia --project=. scripts/generate_pmcmc_seitl.jl
#
# Takes around 15 hours at 256 particles, so run it on the cluster
# (scripts/run_pmcmc.slurm) rather than on a laptop.

include(joinpath(@__DIR__, "pmmh_setup.jl"))

Random.seed!(5678)  # Different seed from SEIT4L

chain = run_pmmh(pmmh_seitl(flu_observations(), N_PARTICLES), "SEITL")

output_path = datadir("pmcmc_seitl_chain.csv")
println("Saving SEITL chain to $output_path")
save_chain_csv(chain, output_path)

print_diagnostics(chain, "SEITL")

println("\n" * "=" ^ 60)
println("Done!")
println("=" ^ 60)
