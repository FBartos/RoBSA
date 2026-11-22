context("(1) JAGS module functionality")
skip_on_cran()

test_that("JAGS 5 loads the module and samples the exponential survival model", {

  RoBSA:::.load_RoBSA_module()
  expect_true(rjags::jags.version() >= "5.0.0")
  expect_true("RoBSA" %in% rjags::list.modules())

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

  connection <- textConnection(model_syntax)
  on.exit(close(connection))
  model <- rjags::jags.model(file = connection, data = data, quiet = TRUE)
  fit   <- rjags::jags.samples(model = model, variable.names = "mu", n.iter = 100, quiet = TRUE, progress.bar = "none")



  expect_equal(summary(fit)[1], "100")
  expect_equal(summary(fit)[2], "mcarray")
  expect_equal(summary(fit)[3], "numeric")
  expect_equal(unname(dim(fit$mu)), c(1, 100, 1))
})
