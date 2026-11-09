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

  expect_equal(
    capture_output_lines(fit, print = TRUE, width = 150),
    c("Call:"                                                                ,
      "RoBSA(formula = Surv(time, cens) ~ x_c + x_d, data = df, chains = 1, ",
      "    sample = 100, burnin = 50, autofit = FALSE, seed = 1)"            ,
      ""                                                                     ,
      "Estimates:"                                                           ,
      "       x_c     x_d[B] "                                               ,
      "-0.1092297 -0.1960542 "                                               ))

  expect_equal(
    capture_output_lines(summary(fit), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                              ,
      "RoBSA(formula = Surv(time, cens) ~ x_c + x_d, data = df, chains = 1, "                                                                                              ,
      "    sample = 100, burnin = 50, autofit = FALSE, seed = 1)"                                                                                                          ,
      ""                                                                                                                                                                   ,
      "Robust Bayesian survival analysis"                                                                                                                                  ,
      "Distributions summary:"                                                                                                                                             ,
      "            Models Prior prob. Post. prob. Inclusion BF"                                                                                                            ,
      "exp-aft       4/20       0.200       0.361        2.263"                                                                                                            ,
      "weibull-aft   4/20       0.200       0.158        0.752"                                                                                                            ,
      "lnorm-aft     4/20       0.200       0.141        0.655"                                                                                                            ,
      "llogis-aft    4/20       0.200       0.153        0.721"                                                                                                            ,
      "gamma-aft     4/20       0.200       0.187        0.920"                                                                                                            ,
      ""                                                                                                                                                                   ,
      "Components summary:"                                                                                                                                                ,
      "    Models Prior prob. Post. prob. Inclusion BF"                                                                                                                    ,
      "x_c  10/20       0.500       0.478        0.914"                                                                                                                    ,
      "x_d  10/20       0.500       0.450        0.819"                                                                                                                    ,
      ""                                                                                                                                                                   ,
      "Model-averaged estimates:"                                                                                                                                          ,
      "         Mean Median  0.025 0.975"                                                                                                                                  ,
      "x_c    -0.109  0.000 -1.482 0.995"                                                                                                                                  ,
      "x_d[B] -0.196  0.000 -1.612 0.664"                                                                                                                                  ,
      "\033[0;31mThe continuous predictors 'x_c' are not scaled. Note that extra care need to be taken when specifying prior distributions for unscaled predictors.\033[0m",
      "\033[0;31mModel (1): ESS 64 is lower than the set target (500). Failing parameters: 'mu_intercept'.\033[0m"                                                         ,
      "\033[0;31mModel (2): ESS 42 is lower than the set target (500). Failing parameters: 'mu_intercept', 'aux'.\033[0m"                                                  ,
      "\033[0;31mModel (3): ESS 43 is lower than the set target (500). Failing parameters: 'mu_intercept', 'aux'.\033[0m"                                                  ,
      "\033[0;31mModel (4): ESS 59 is lower than the set target (500). Failing parameters: 'mu_intercept', 'aux'.\033[0m"                                                  ,
      "\033[0;31mThere were another 15 warnings. To see all warnings call 'check_RoBSA(fit)'.\033[0m"))

})

test_that("BayesTools version check accepts development versions", {

  expect_true(RoBSA:::.check_BayesTools(package_version("1.1.0"),      package_version("0.4.0")))
  expect_true(RoBSA:::.check_BayesTools(package_version("1.1.0.9001"), package_version("0.4.0.9001")))
  expect_true(RoBSA:::.check_BayesTools(package_version("1.1.0"),      package_version("0.10.0")))

  expect_warning(
    ok <- RoBSA:::.check_BayesTools(package_version("1.1.0.9001"), package_version("0.3.0")),
    "RoBSA version 1.1.0.9001 requires BayesTools version higher or equal 0.4.0"
  )
  expect_false(ok)
})

test_that("Parallel fitting uses at least one worker", {

  expect_equal(RoBSA:::.parallel_workers(7, 3), 2)
  expect_equal(RoBSA:::.parallel_workers(3, 3), 1)
  expect_equal(RoBSA:::.parallel_workers(2, 3), 1)
  expect_equal(RoBSA:::.parallel_workers(0, 3), 1)
})
