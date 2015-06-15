## library('fitR')

fitR.dir <- switch(Sys.info()[["user"]],
	Tonton = "~/edu/Fit_course/fitR",## Anton
	seb ="~/teaching/fitR"
	) 

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-fitmodel.r")))
save(SIR, file = path.expand(paste0(fitR.dir, "/", "data/SIR.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-fitmodel-stochastic.r")))
save(SIR_stoch, file = path.expand(paste0(fitR.dir, "/", "data/SIR_stoch.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-fitmodel-log.r")))
save(SIR_exp, file = path.expand(paste0(fitR.dir, "/", "data/SIR_exp.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-fitmodel-reporting.r")))
save(SIR_reporting, file = path.expand(paste0(fitR.dir, "/", "data/SIR_reporting.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-SEITL-deter.r")))
save(SEITL_deter, file = path.expand(paste0(fitR.dir, "/", "data/SEITL_deter.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-SEITL-sto.r")))
save(SEITL_stoch, file = path.expand(paste0(fitR.dir, "/", "data/SEITL_stoch.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-SEIT2L-deter.r")))
save(SEIT2L_deter, file = path.expand(paste0(fitR.dir, "/", "data/SEIT2L_deter.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-SEIT2L-sto.r")))
save(SEIT2L_stoch, file = path.expand(paste0(fitR.dir, "/", "data/SEIT2L_stoch.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-SEIT4L-deter.r")))
save(SEIT4L_deter, file = path.expand(paste0(fitR.dir, "/", "data/SEIT4L_deter.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-SEIT4L-sto.r")))
save(SEIT4L_stoch, file = path.expand(paste0(fitR.dir, "/", "data/SEIT4L_stoch.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-SEITL-pomp.r")))
save(SEITL_pomp, file = path.expand(paste0(fitR.dir, "/", "data/SEITL_pomp.rdata")))

source(path.expand(paste0(fitR.dir, "/", "inst/examples/example-SEIT2L-pomp.r")))
save(SEIT2L_pomp, file = path.expand(paste0(fitR.dir, "/", "data/SEIT2L_pomp.rdata")))


