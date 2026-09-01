# Generate both pre-computed PMMH chains loaded by sessions/pmcmc.qmd.
# Run with: julia --project=. scripts/generate_pmcmc_chains.jl
#
# The two models take around 15 hours each at 256 particles, which is why their
# chains are committed to data/ and the session loads those. Run generate_pmcmc_seitl.jl
# and generate_pmcmc_seit4l.jl side by side to halve the wall clock if regenerating.

include(joinpath(@__DIR__, "pmmh_setup.jl"))

obs = flu_observations()

Random.seed!(1234)
chain_seit4l = run_pmmh(pmmh_seit4l(obs, N_PARTICLES), "SEIT4L")

output_path = datadir("pmcmc_seit4l_chain.csv")
println("Saving SEIT4L chain to $output_path")
save_chain_csv(chain_seit4l, output_path)

print_diagnostics(chain_seit4l, "SEIT4L")

Random.seed!(5678)  # Different seed from SEIT4L
chain_seitl = run_pmmh(pmmh_seitl(obs, N_PARTICLES), "SEITL")

output_path_seitl = datadir("pmcmc_seitl_chain.csv")
println("Saving SEITL chain to $output_path_seitl")
save_chain_csv(chain_seitl, output_path_seitl)

print_diagnostics(chain_seitl, "SEITL")

println("\n" * "=" ^ 60)
println("All done!")
println("=" ^ 60)
