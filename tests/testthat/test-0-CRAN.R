context("(0) Basic tests for CRAN")
# These are just a very rudimentary tests that don't require time or saved files.
# The full range of tests is run locally.


test_that("Basic functionality works", {

  RoBSA:::.load_RoBSA_module()

  # fit a default model
  df <- data.frame(
    time = c(0.5, 1.0, 0.2, 0.8, 0.1, 0.1, 0.8, 0.4, 0.2, 0.8),
    cens = c(  0,   1,   1,   0,   0,   1,   0,   0,   1,   1),
    x_c  = c(0.2, 0.8, 0.1, 0.3, 0.5, 0.4, 0.2, 0.1, 0.8, 0.9),
    x_d  = c("A", "A", "A", "B", "B", "B", "A", "A", "B", "B")
  )

  fit <- suppressWarnings(RoBSA(Surv(time, cens) ~ x_c + x_d, data = df, chains = 1, burnin = 50, sample = 100, autofit = FALSE, seed = 1))

  expect_equal(TRUE, is.RoBSA(fit))
  for(model in fit$models){
    expect_s3_class(model$fit, "BayesTools_fit")
    expect_equal(coda::niter(coda::as.mcmc.list(model$fit)), 100)
  }

  expect_equal(
    encodeString(capture_output_lines(fit, print = TRUE, width = 150), quote = '"'),
    readLines(test_path("reference", "basic-print.txt"))
  )

  expect_equal(
    encodeString(capture_output_lines(summary(fit), print = TRUE, width = 150), quote = '"'),
    readLines(test_path("reference", "basic-summary.txt"))
  )

})

test_that("BayesTools version check accepts development versions", {

  expect_true(RoBSA:::.check_BayesTools(package_version("1.1.0"),      package_version("0.4.0")))
  expect_true(RoBSA:::.check_BayesTools(package_version("1.1.0.9001"), package_version("0.4.0.9001")))
  expect_true(RoBSA:::.check_BayesTools(package_version("1.1.0"),      package_version("0.10.0")))
  expect_true(RoBSA:::.check_BayesTools(package_version("1.1.0.9003"), package_version("0.4.0.9013")))

  expect_warning(
    ok <- RoBSA:::.check_BayesTools(package_version("1.1.0.9001"), package_version("0.3.0")),
    "RoBSA version 1.1.0.9001 requires BayesTools version higher or equal 0.4.0"
  )
  expect_false(ok)
  expect_warning(
    ok <- RoBSA:::.check_BayesTools(package_version("1.1.0.9003"), package_version("0.4.0.9012")),
    "0.4.0.9013", fixed = TRUE
  )
  expect_false(ok)
})

test_that("Parallel fitting uses at least one worker", {

  expect_equal(RoBSA:::.parallel_workers(7, 3), 2)
  expect_equal(RoBSA:::.parallel_workers(3, 3), 1)
  expect_equal(RoBSA:::.parallel_workers(2, 3), 1)
  expect_equal(RoBSA:::.parallel_workers(0, 3), 1)
})
