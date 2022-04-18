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
      "-0.1708049 -0.1788191 "))

  expect_equal(
    capture_output_lines(summary(fit), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                              ,
      "RoBSA(formula = Surv(time, cens) ~ x_c + x_d, data = df, chains = 1, "                                                                                              ,
      "    sample = 100, burnin = 50, autofit = FALSE, seed = 1)"                                                                                                          ,
      ""                                                                                                                                                                   ,
      "Robust Bayesian survival analysis"                                                                                                                                  ,
      "Distributions summary:"                                                                                                                                             ,
      "            Models Prior prob. Post. prob. Inclusion BF"                                                                                                            ,
      "exp-aft       4/20       0.200       0.361        2.257"                                                                                                            ,
      "weibull-aft   4/20       0.200       0.156        0.741"                                                                                                            ,
      "lnorm-aft     4/20       0.200       0.146        0.684"                                                                                                            ,
      "llogis-aft    4/20       0.200       0.146        0.682"                                                                                                            ,
      "gamma-aft     4/20       0.200       0.191        0.946"                                                                                                            ,
      ""                                                                                                                                                                   ,
      "Components summary:"                                                                                                                                                ,
      "         Models Prior prob. Post. prob. Inclusion BF"                                                                                                               ,
      "(mu) x_c  10/20       0.500       0.473        0.898"                                                                                                               ,
      "(mu) x_d  10/20       0.500       0.430        0.756"                                                                                                               ,
      ""                                                                                                                                                                   ,
      "Model-averaged estimates:"                                                                                                                                          ,
      "         Mean Median  0.025 0.975"                                                                                                                                  ,
      "x_c    -0.171  0.000 -1.674 0.868"                                                                                                                                  ,
      "x_d[B] -0.179  0.000 -1.432 0.620"                                                                                                                                  ,
      "\033[0;31mThe continuous predictors 'x_c' are not scaled. Note that extra care need to be taken when specifying prior distributions for unscaled predictors.\033[0m",
      "\033[0;31mModel (1): ESS 84 is lower than the set target (500).\033[0m"                                                                                             ,
      "\033[0;31mModel (2): ESS 54 is lower than the set target (500).\033[0m"                                                                                             ,
      "\033[0;31mModel (3): ESS 17 is lower than the set target (500).\033[0m"                                                                                             ,
      "\033[0;31mModel (4): ESS 34 is lower than the set target (500).\033[0m"                                                                                             ,
      "\033[0;31mThere were another 15 warnings. To see all warnings call 'check_RoBSA(fit)'.\033[0m"))

})
