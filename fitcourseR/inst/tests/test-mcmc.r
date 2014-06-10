context("mcmc")

test_that("mcmcMH for deterministic SEIT2L model",{

	SEIT2L <- SEIT2L_createModelTdC(deterministic=TRUE, verbose=FALSE) 

	theta.init <- SEIT2L$theta

	suppressMessages(ans <- mcmcMH(target=targetPosterior, target.args=list(log.prior=SEIT2L$log.prior, marginal.log.likelihood= marginalLogLikelihoodDeterministic, marginal.log.likelihood.args=list(fitmodel=SEIT2L)), theta.init=theta.init, gaussian.proposal=SEIT2L$gaussian.proposal, n.iterations=100, adapt.size.start=10, adapt.size.cooling=0.99, adapt.shape.start=10, print.info.every=NULL))

	expect_true(is.data.frame(ans$trace))
	expect_true(is.numeric(ans$acceptance.rate))
	expect_true(is.matrix(ans$covmat.empirical))

})

test_that("mcmcMH for stochastic SEIT2L model",{

	SEIT2L <- SEIT2L_createModelTdC(deterministic=FALSE, verbose=FALSE) 

	theta.init <- SEIT2L$theta

	suppressMessages(ans <- mcmcMH(target=targetPosterior, target.args=list(log.prior=SEIT2L$log.prior, marginal.log.likelihood= marginalLogLikelihoodStochastic, marginal.log.likelihood.args=list(fitmodel=SEIT2L,n.particles=20,n.cores=1)), theta.init=theta.init, gaussian.proposal=SEIT2L$gaussian.proposal, n.iterations=2, adapt.size.start=10, adapt.size.cooling=0.99, adapt.shape.start=10, print.info.every=NULL))

	expect_true(is.data.frame(ans$trace))
	expect_true(is.numeric(ans$acceptance.rate))
	expect_true(is.matrix(ans$covmat.empirical))

})


test_that("mcmcMH ABC for deterministic SEIT2L model",{

	SEIT2L <- SEIT2L_createModelTdC(deterministic=TRUE, verbose=FALSE) 

	theta.init <- SEIT2L$theta

	suppressMessages(ans <- mcmcMH(target=targetPosteriorABC, target.args=list(fitmodel=SEIT2L,epsilon=1), theta.init=theta.init, gaussian.proposal=SEIT2L$gaussian.proposal, n.iterations=100, adapt.size.start=10, adapt.size.cooling=0.99, adapt.shape.start=10, print.info.every=NULL))

	expect_true(is.data.frame(ans$trace))
	expect_true(is.numeric(ans$acceptance.rate))
	expect_true(is.matrix(ans$covmat.empirical))

})
