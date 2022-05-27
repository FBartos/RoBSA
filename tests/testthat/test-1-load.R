context("(1) JAGS module functionality")
skip_on_cran()

### the regression parameterized distributions can be loaded and run
test_that("Module can be loaded and the one-sided normal distribution works", {

  module_location <- gsub('/$','', file.path(.libPaths(), "RoBSA", 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
  sapply(module_location, function(path) rjags::load.module("RoBSA", path = path))

  model_syntax <- "
  model{
    mu ~ dnorm(0, 1)
    for(i in 1:K){
      t[i] ~ exp_aft_event(mu)
    }
  }"

  data <- list(
    t      = rgamma(10, 1, 1),
    K      = 10
  )

  RoBSA:::.load_RoBSA_module()
  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
  fit   <- rjags::jags.samples(model = model, variable.names = "mu", n.iter = 100, quiet = TRUE, progress.bar = "none")



  expect_equal(summary(fit)[1], "100")
  expect_equal(summary(fit)[2], "mcarray")
  expect_equal(summary(fit)[3], "numeric")
  expect_equal(unname(dim(fit$mu)), c(1, 100, 1))
})
