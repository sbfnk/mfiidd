context("classes")

test_that("constructor for SEITL model",{

	SEITL_det <- SEITL_createModelTdC(deterministic=TRUE,verbose=FALSE)
	expect_true(inherits(SEITL_det,"fitmodel"))


	SEITL_sto <- SEITL_createModelTdC(deterministic=FALSE,verbose=FALSE)
	expect_true(inherits(SEITL_sto,"fitmodel"))


})