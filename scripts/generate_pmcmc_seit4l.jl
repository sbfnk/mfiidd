# Generate the pre-computed SEIT4L PMMH chain loaded by sessions/pmcmc.qmd.
# Run with: julia --project=. scripts/generate_pmcmc_seit4l.jl
#
# Takes around 14 hours at 256 particles, which is why the chains it produces
# are committed to data/ and the session loads those instead of running this.

include(joinpath(@__DIR__, "pmmh_setup.jl"))

Random.seed!(1234)

chain = run_pmmh(pmmh_seit4l(flu_observations(), N_PARTICLES), "SEIT4L")

output_path = datadir("pmcmc_seit4l_chain.csv")
println("Saving SEIT4L chain to $output_path")
save_chain_csv(chain, output_path)

print_diagnostics(chain, "SEIT4L")

println("\n" * "=" ^ 60)
println("Done!")
println("=" ^ 60)
