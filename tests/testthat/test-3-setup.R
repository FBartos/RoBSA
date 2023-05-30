context("(3) Model setup")
skip_on_cran()

# test model preview
test_that("Model preview & set-up works", {

  df <- data.frame(
    time   = rgamma(60, 1, 1),
    event  = rbinom(60, 1, 0.5),
    x_cont = rnorm(60),
    x_bin  = as.factor(rbinom(60, 1, .5)),
    x_fac3 = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac4 = factor(rep(c("A", "B", "C", "D"), 15), levels = c("A", "B", "C", "D"))
  )


  ### intercept only
  expect_equal(
    capture_output_lines(check_setup(Surv(time = time, event = event) ~ 1, data = df, models = FALSE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)",
      "Distributions summary:"                    ,
      "            Models Prior prob."            ,
      "exp-aft        1/5       0.200"            ,
      "weibull-aft    1/5       0.200"            ,
      "lnorm-aft      1/5       0.200"            ,
      "llogis-aft     1/5       0.200"            ,
      "gamma-aft      1/5       0.200"
    )
  )

  expect_equal(
    capture_output_lines(check_setup(Surv(time = time, event = event) ~ 1, data = df, models = TRUE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)"                          ,
      "Models overview:"                                                    ,
      " Model Distribution Prior Intercept    Prior Auxiliary   Prior prob.",
      "     1      exp-aft    Normal(0, 5)                 None       0.200",
      "     2  weibull-aft    Normal(0, 5) Normal(0, 1)[0, Inf]       0.200",
      "     3    lnorm-aft    Normal(0, 5) Normal(0, 1)[0, Inf]       0.200",
      "     4   llogis-aft    Normal(0, 5) Normal(0, 1)[0, Inf]       0.200",
      "     5    gamma-aft    Normal(0, 5) Normal(0, 1)[0, Inf]       0.200"
    )
  )

  ### test all (default)
  expect_equal(
    capture_output_lines(check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                                     data = df, models = FALSE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)",
      "Distributions summary:"                    ,
      "            Models Prior prob."            ,
      "exp-aft       4/20       0.200"            ,
      "weibull-aft   4/20       0.200"            ,
      "lnorm-aft     4/20       0.200"            ,
      "llogis-aft    4/20       0.200"            ,
      "gamma-aft     4/20       0.200"            ,
      ""                                          ,
      "Components summary:"                       ,
      "       Models Prior prob."                 ,
      "x_cont  10/20       0.500"                 ,
      "x_bin   10/20       0.500"
    )
  )

  expect_equal(
    capture_output_lines(check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                                     data = df, models = TRUE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)",
      "Models overview:"                                                                                                  ,
      " Model Distribution Prior Intercept    Prior Auxiliary   Prior x_cont            Prior x_bin           Prior prob.",
      "     1      exp-aft    Normal(0, 5)                 None     Spike(0)     treatment contrast: Spike(0)       0.050",
      "     2  weibull-aft    Normal(0, 5) Normal(0, 1)[0, Inf]     Spike(0)     treatment contrast: Spike(0)       0.050",
      "     3    lnorm-aft    Normal(0, 5) Normal(0, 1)[0, Inf]     Spike(0)     treatment contrast: Spike(0)       0.050",
      "     4   llogis-aft    Normal(0, 5) Normal(0, 1)[0, Inf]     Spike(0)     treatment contrast: Spike(0)       0.050",
      "     5    gamma-aft    Normal(0, 5) Normal(0, 1)[0, Inf]     Spike(0)     treatment contrast: Spike(0)       0.050",
      "     6      exp-aft    Normal(0, 5)                 None Normal(0, 1)     treatment contrast: Spike(0)       0.050",
      "     7  weibull-aft    Normal(0, 5) Normal(0, 1)[0, Inf] Normal(0, 1)     treatment contrast: Spike(0)       0.050",
      "     8    lnorm-aft    Normal(0, 5) Normal(0, 1)[0, Inf] Normal(0, 1)     treatment contrast: Spike(0)       0.050",
      "     9   llogis-aft    Normal(0, 5) Normal(0, 1)[0, Inf] Normal(0, 1)     treatment contrast: Spike(0)       0.050",
      "    10    gamma-aft    Normal(0, 5) Normal(0, 1)[0, Inf] Normal(0, 1)     treatment contrast: Spike(0)       0.050",
      "    11      exp-aft    Normal(0, 5)                 None     Spike(0) treatment contrast: Normal(0, 1)       0.050",
      "    12  weibull-aft    Normal(0, 5) Normal(0, 1)[0, Inf]     Spike(0) treatment contrast: Normal(0, 1)       0.050",
      "    13    lnorm-aft    Normal(0, 5) Normal(0, 1)[0, Inf]     Spike(0) treatment contrast: Normal(0, 1)       0.050",
      "    14   llogis-aft    Normal(0, 5) Normal(0, 1)[0, Inf]     Spike(0) treatment contrast: Normal(0, 1)       0.050",
      "    15    gamma-aft    Normal(0, 5) Normal(0, 1)[0, Inf]     Spike(0) treatment contrast: Normal(0, 1)       0.050",
      "    16      exp-aft    Normal(0, 5)                 None Normal(0, 1) treatment contrast: Normal(0, 1)       0.050",
      "    17  weibull-aft    Normal(0, 5) Normal(0, 1)[0, Inf] Normal(0, 1) treatment contrast: Normal(0, 1)       0.050",
      "    18    lnorm-aft    Normal(0, 5) Normal(0, 1)[0, Inf] Normal(0, 1) treatment contrast: Normal(0, 1)       0.050",
      "    19   llogis-aft    Normal(0, 5) Normal(0, 1)[0, Inf] Normal(0, 1) treatment contrast: Normal(0, 1)       0.050",
      "    20    gamma-aft    Normal(0, 5) Normal(0, 1)[0, Inf] Normal(0, 1) treatment contrast: Normal(0, 1)       0.050"
    )
  )

  ### test one (custom priors)
  expect_equal(
    capture_output_lines(check_setup(Surv(time = time, event = event) ~ x_cont + x_bin, test_predictors = "x_bin",
                                     priors = list(
                                       "x_bin"  = prior_factor("mnormal", list(0, .33), contrast = "orthonormal"),
                                       "x_cont" = prior("normal", list(0, .25))),
                                     data = df, models = FALSE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)",
      "Distributions summary:"                    ,
      "            Models Prior prob."            ,
      "exp-aft       2/10       0.200"            ,
      "weibull-aft   2/10       0.200"            ,
      "lnorm-aft     2/10       0.200"            ,
      "llogis-aft    2/10       0.200"            ,
      "gamma-aft     2/10       0.200"            ,
      ""                                          ,
      "Components summary:"                       ,
      "      Models Prior prob."                  ,
      "x_bin   5/10       0.500"
    )
  )

  expect_equal(
    capture_output_lines(check_setup(Surv(time = time, event = event) ~ x_cont + x_bin, test_predictors = "x_bin",
                                     priors = list(
                                       "x_bin"  = prior_factor("mnormal", list(0, .33), contrast = "orthonormal"),
                                       "x_cont" = prior("normal", list(0, .25))),
                                     data = df, models = TRUE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)"                                                                                   ,
      "Models overview:"                                                                                                             ,
      " Model Distribution Prior Intercept    Prior Auxiliary     Prior x_cont                 Prior x_bin               Prior prob.",
      "     1      exp-aft    Normal(0, 5)                 None  Normal(0, 0.25)         orthonormal contrast: mSpike(0)       0.100",
      "     2  weibull-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)         orthonormal contrast: mSpike(0)       0.100",
      "     3    lnorm-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)         orthonormal contrast: mSpike(0)       0.100",
      "     4   llogis-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)         orthonormal contrast: mSpike(0)       0.100",
      "     5    gamma-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)         orthonormal contrast: mSpike(0)       0.100",
      "     6      exp-aft    Normal(0, 5)                 None  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)       0.100",
      "     7  weibull-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)       0.100",
      "     8    lnorm-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)       0.100",
      "     9   llogis-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)       0.100",
      "    10    gamma-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)       0.100"
    )
  )

  ### test all ( + change default priors)
  expect_equal(
    capture_output_lines( check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                                      distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
                                      distributions_weights = c(3, 1, 1, 1),
                                      prior_beta_null   = prior("uniform", list(-0.1, 0.1)),
                                      prior_factor_null = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                                      prior_beta        = prior("gamma", list(10, 10)),
                                      prior_factor      = prior_factor("gamma", list(20, 20), contrast = "treatment"),
                                      prior_intercept = list(
                                        "gamma-aft"   = prior("normal", list(1, 1)),
                                        "weibull-aft" = prior("normal", list(2, 1)),
                                        "lnorm-aft"   = prior("normal", list(1, 3)),
                                        "llogis-aft"  = prior("normal", list(1, 4))
                                      ),
                                      data = df, models = FALSE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)",
      "Distributions summary:"                    ,
      "            Models Prior prob."            ,
      "gamma-aft     4/16       0.500"            ,
      "weibull-aft   4/16       0.167"            ,
      "lnorm-aft     4/16       0.167"            ,
      "llogis-aft    4/16       0.167"            ,
      ""                                          ,
      "Components summary:"                       ,
      "       Models Prior prob."                 ,
      "x_cont   8/16       0.500"                 ,
      "x_bin    8/16       0.500"
    )
  )

  expect_equal(
    capture_output_lines( check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                                      distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
                                      distributions_weights = c(3, 1, 1, 1),
                                      prior_beta_null   = prior("uniform", list(-0.1, 0.1)),
                                      prior_factor_null = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                                      prior_beta        = prior("gamma", list(10, 10)),
                                      prior_factor      = prior_factor("gamma", list(20, 20), contrast = "treatment"),
                                      prior_intercept = list(
                                        "gamma-aft"   = prior("normal", list(1, 1)),
                                        "weibull-aft" = prior("normal", list(2, 1)),
                                        "lnorm-aft"   = prior("normal", list(1, 3)),
                                        "llogis-aft"  = prior("normal", list(1, 4))
                                      ),
                                      data = df, models = TRUE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)"                                                                                     ,
      "Models overview:"                                                                                                               ,
      " Model Distribution Prior Intercept    Prior Auxiliary      Prior x_cont                  Prior x_bin               Prior prob.",
      "     1    gamma-aft    Normal(1, 1) Normal(0, 1)[0, Inf] Uniform(-0.1, 0.1)  treatment contrast: Uniform(-0.2, 0.2)       0.125",
      "     2  weibull-aft    Normal(2, 1) Normal(0, 1)[0, Inf] Uniform(-0.1, 0.1)  treatment contrast: Uniform(-0.2, 0.2)       0.042",
      "     3    lnorm-aft    Normal(1, 3) Normal(0, 1)[0, Inf] Uniform(-0.1, 0.1)  treatment contrast: Uniform(-0.2, 0.2)       0.042",
      "     4   llogis-aft    Normal(1, 4) Normal(0, 1)[0, Inf] Uniform(-0.1, 0.1)  treatment contrast: Uniform(-0.2, 0.2)       0.042",
      "     5    gamma-aft    Normal(1, 1) Normal(0, 1)[0, Inf]      Gamma(10, 10)  treatment contrast: Uniform(-0.2, 0.2)       0.125",
      "     6  weibull-aft    Normal(2, 1) Normal(0, 1)[0, Inf]      Gamma(10, 10)  treatment contrast: Uniform(-0.2, 0.2)       0.042",
      "     7    lnorm-aft    Normal(1, 3) Normal(0, 1)[0, Inf]      Gamma(10, 10)  treatment contrast: Uniform(-0.2, 0.2)       0.042",
      "     8   llogis-aft    Normal(1, 4) Normal(0, 1)[0, Inf]      Gamma(10, 10)  treatment contrast: Uniform(-0.2, 0.2)       0.042",
      "     9    gamma-aft    Normal(1, 1) Normal(0, 1)[0, Inf] Uniform(-0.1, 0.1)       treatment contrast: Gamma(20, 20)       0.125",
      "    10  weibull-aft    Normal(2, 1) Normal(0, 1)[0, Inf] Uniform(-0.1, 0.1)       treatment contrast: Gamma(20, 20)       0.042",
      "    11    lnorm-aft    Normal(1, 3) Normal(0, 1)[0, Inf] Uniform(-0.1, 0.1)       treatment contrast: Gamma(20, 20)       0.042",
      "    12   llogis-aft    Normal(1, 4) Normal(0, 1)[0, Inf] Uniform(-0.1, 0.1)       treatment contrast: Gamma(20, 20)       0.042",
      "    13    gamma-aft    Normal(1, 1) Normal(0, 1)[0, Inf]      Gamma(10, 10)       treatment contrast: Gamma(20, 20)       0.125",
      "    14  weibull-aft    Normal(2, 1) Normal(0, 1)[0, Inf]      Gamma(10, 10)       treatment contrast: Gamma(20, 20)       0.042",
      "    15    lnorm-aft    Normal(1, 3) Normal(0, 1)[0, Inf]      Gamma(10, 10)       treatment contrast: Gamma(20, 20)       0.042",
      "    16   llogis-aft    Normal(1, 4) Normal(0, 1)[0, Inf]      Gamma(10, 10)       treatment contrast: Gamma(20, 20)       0.042"
    )
  )

  ### test none + custom priors + interaction
  expect_equal(
    capture_output_lines(check_setup(Surv(time = time, event = event) ~ x_cont * x_bin, test_predictors = FALSE,
                                     priors = list(
                                       "x_bin"        = prior_factor("mnormal", list(0, .33), contrast = "orthonormal"),
                                       "x_cont"       = prior("normal", list(0, .25)),
                                       "x_cont:x_bin" = prior_factor("mnormal", list(0, .10), contrast = "orthonormal")),
                                     data = df, models = FALSE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)",
      "Distributions summary:"                    ,
      "            Models Prior prob."            ,
      "exp-aft        1/5       0.200"            ,
      "weibull-aft    1/5       0.200"            ,
      "lnorm-aft      1/5       0.200"            ,
      "llogis-aft     1/5       0.200"            ,
      "gamma-aft      1/5       0.200"
    )
  )

  expect_equal(
    capture_output_lines(check_setup(Surv(time = time, event = event) ~ x_cont * x_bin, test_predictors = FALSE,
                                     priors = list(
                                       "x_bin"        = prior_factor("mnormal", list(0, .33), contrast = "orthonormal"),
                                       "x_cont"       = prior("normal", list(0, .25)),
                                       "x_cont:x_bin" = prior_factor("mnormal", list(0, .10), contrast = "orthonormal")),
                                     data = df, models = TRUE, rescale_data = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian survival analysis (set-up)"                                                                       ,
      "Models overview:"                                                                                                 ,
      " Model Distribution Prior Intercept    Prior Auxiliary     Prior x_cont                 Prior x_bin              ",
      "     1      exp-aft    Normal(0, 5)                 None  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)",
      "     2  weibull-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)",
      "     3    lnorm-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)",
      "     4   llogis-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)",
      "     5    gamma-aft    Normal(0, 5) Normal(0, 1)[0, Inf]  Normal(0, 0.25)  orthonormal contrast: mNormal(0, 0.33)",
      "           Prior x_cont:x_bin           Prior prob."                                                              ,
      "  orthonormal contrast: mNormal(0, 0.1)       0.200"                                                              ,
      "  orthonormal contrast: mNormal(0, 0.1)       0.200"                                                              ,
      "  orthonormal contrast: mNormal(0, 0.1)       0.200"                                                              ,
      "  orthonormal contrast: mNormal(0, 0.1)       0.200"                                                              ,
      "  orthonormal contrast: mNormal(0, 0.1)       0.200"
    )
  )


})

test_that("Prior input checks work", {

  df <- data.frame(
    time   = rgamma(60, 1, 1),
    event  = rbinom(60, 1, 0.5),
    x_cont = rnorm(60),
    x_bin  = as.factor(rbinom(60, 1, .5)),
    x_fac3 = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac4 = factor(rep(c("A", "B", "C", "D"), 15), levels = c("A", "B", "C", "D"))
  )


  expect_error(check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                             distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
                             distributions_weights = c(3, 1, 1, 1),
                             prior_beta_null   = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                             prior_factor_null = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                             prior_beta        = prior("gamma", list(10, 10)),
                             prior_factor      = prior_factor("gamma", list(20, 20), contrast = "treatment"),
                             prior_intercept = list(
                               "gamma-aft"   = prior("normal", list(1, 1)),
                               "weibull-aft" = prior("normal", list(2, 1)),
                               "lnorm-aft"   = prior("normal", list(1, 3)),
                               "llogis-aft"  = prior("normal", list(1, 4))
                             ),
                             data = df, models = FALSE, rescale_data = TRUE), "The default prior for predictors (null) is not a valid prior distribution.", fixed = TRUE)

  expect_error(check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                             distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
                             distributions_weights = c(3, 1, 1, 1),
                             prior_beta_null   = prior("spike", list(0)),
                             prior_factor_null = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                             prior_beta        = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                             prior_factor      = prior_factor("gamma", list(20, 20), contrast = "treatment"),
                             prior_intercept = list(
                               "gamma-aft"   = prior("normal", list(1, 1)),
                               "weibull-aft" = prior("normal", list(2, 1)),
                               "lnorm-aft"   = prior("normal", list(1, 3)),
                               "llogis-aft"  = prior("normal", list(1, 4))
                             ),
                             data = df, models = FALSE, rescale_data = TRUE), "The default prior for predictors (alt) is not a valid prior distribution.", fixed = TRUE)

  expect_error(check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                             distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
                             distributions_weights = c(3, 1, 1, 1),
                             prior_beta_null   = prior("uniform", list(-0.1, 0.1)),
                             prior_factor_null = prior("uniform", list(-0.2, 0.2)),
                             prior_beta        = prior("gamma", list(10, 10)),
                             prior_factor      = prior_factor("gamma", list(20, 20), contrast = "treatment"),
                             prior_intercept = list(
                               "gamma-aft"   = prior("normal", list(1, 1)),
                               "weibull-aft" = prior("normal", list(2, 1)),
                               "lnorm-aft"   = prior("normal", list(1, 3)),
                               "llogis-aft"  = prior("normal", list(1, 4))
                             ),
                             data = df, models = FALSE, rescale_data = TRUE), "The default prior for factors (null) is not a valid prior distribution.", fixed = TRUE)

  expect_error(  check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                             distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
                             distributions_weights = c(3, 1, 1, 1),
                             prior_beta_null   = prior("uniform", list(-0.1, 0.1)),
                             prior_factor_null = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                             prior_beta        = prior("gamma", list(10, 10)),
                             prior_factor      = prior("gamma", list(10, 10)),
                             prior_intercept = list(
                               "gamma-aft"   = prior("normal", list(1, 1)),
                               "weibull-aft" = prior("normal", list(2, 1)),
                               "lnorm-aft"   = prior("normal", list(1, 3)),
                               "llogis-aft"  = prior("normal", list(1, 4))
                             ),
                             data = df, models = FALSE, rescale_data = TRUE), "The default prior for factors (alt) is not a valid prior distribution.", fixed = TRUE)


  expect_error(  check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                             distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
                             distributions_weights = c(3, 1, 1, 1),
                             prior_beta_null   = prior("uniform", list(-0.1, 0.1)),
                             prior_factor_null = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                             prior_beta        = prior("gamma", list(10, 10)),
                             prior_factor      = prior_factor("gamma", list(20, 20), contrast = "treatment"),
                             prior_intercept = list(
                               "gamma-aft"   = prior_factor("gamma", list(20, 20), contrast = "treatment"),
                               "weibull-aft" = prior("normal", list(2, 1)),
                               "lnorm-aft"   = prior("normal", list(1, 3)),
                               "llogis-aft"  = prior("normal", list(1, 4))
                             ),
                             data = df, models = FALSE, rescale_data = TRUE), "The default prior for intercepts are not a valid prior distribution.", fixed = TRUE)
  expect_error(  check_setup(Surv(time = time, event = event) ~ x_cont + x_bin,
                             distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
                             distributions_weights = c(3, 1, 1, 1),
                             prior_beta_null   = prior("uniform", list(-0.1, 0.1)),
                             prior_factor_null = prior_factor("uniform", list(-0.2, 0.2), contrast = "treatment"),
                             prior_beta        = prior("gamma", list(10, 10)),
                             prior_factor      = prior_factor("gamma", list(20, 20), contrast = "treatment"),
                             prior_aux = list(
                               "gamma-aft"   = prior("normal", list(1, 1)),
                               "weibull-aft" = prior("normal", list(2, 1)),
                               "lnorm-aft"   = prior("normal", list(1, 3)),
                               "llogis-aft"  = prior_factor("gamma", list(20, 20), contrast = "treatment")
                             ),
                             data = df, models = FALSE, rescale_data = TRUE), "The default prior for auxilary parameters are not a valid prior distribution.", fixed = TRUE)


  expect_error(  RoBSA(Surv(time = time, event = event) ~ x_cont + x_bin, chains = 1, autofit = FALSE, sample = 200, burnin = 100,
                             priors = list(
                               "x_bin"  = prior("normal", list(0, .25)),
                               "x_cont" = prior("normal", list(0, .25))),
                             data = df, rescale_data = TRUE), "Unsupported prior distribution defined for 'x_bin' factor variable. See '?prior_factor' for details.", fixed = TRUE)

  expect_error(  RoBSA(Surv(time = time, event = event) ~ x_cont + x_bin, chains = 1, autofit = FALSE, sample = 200, burnin = 100,
                             priors = list(
                               "x_bin"  = prior_factor("mnormal", list(0, .33), contrast = "orthonormal"),
                               "x_cont" = prior_factor("mnormal", list(0, .33), contrast = "orthonormal")),
                             data = df, rescale_data = TRUE), "Unsupported prior distribution defined for 'x_cont' continuous variable. See '?prior' for details.", fixed = TRUE)

})

test_that("Set autofit control works", {

  expect_error(set_autofit_control(max_Rhat = .99), "Checking 'autofit_control':\n\tThe 'max_Rhat' must be equal or higher than 1.")
  expect_error(set_autofit_control(min_ESS  =  -1), "Checking 'autofit_control':\n\tThe 'min_ESS' must be equal or higher than 0.")
  expect_error(set_autofit_control(max_error=  -1), "Checking 'autofit_control':\n\tThe 'max_error' must be equal or higher than 0.")
  expect_error(set_autofit_control(max_SD_error=  -.1), "Checking 'autofit_control':\n\tThe 'max_SD_error' must be equal or higher than 0.")
  expect_error(set_autofit_control(max_SD_error=  1.1), "Checking 'autofit_control':\n\tThe 'max_SD_error' must be equal or lower than 1.")
  expect_error(set_autofit_control(max_time = list(time = -1, unit = "secs")), "Checking 'autofit_control':\n\tThe 'max_time:time' must be equal or higher than 0.")
  expect_error(set_autofit_control(max_time = list(time = 10, unit = "maps")), "Checking 'autofit_control':\n\tThe 'maps' values are not recognized by the 'max_time:unit' argument.")

  expect_equal(set_autofit_control(), list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(max_Rhat = 1.01),  list(
    max_Rhat      = 1.01,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(min_ESS = 200),  list(
    max_Rhat      = 1.05,
    min_ESS       = 200,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(max_error = 0.01),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = 0.01,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(max_SD_error = 0.01),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = 0.01,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(max_time = list(time = 30, unit = "secs")),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 30, unit = "secs"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(sample_extend = 200),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 200
  ))

})

test_that("Set convergence checks works", {


  expect_error(set_convergence_checks(max_Rhat = .99), "Checking 'convergence_checks':\n\tThe 'max_Rhat' must be equal or higher than 1.")
  expect_error(set_convergence_checks(min_ESS  =  -1), "Checking 'convergence_checks':\n\tThe 'min_ESS' must be equal or higher than 0.")
  expect_error(set_convergence_checks(max_error=  -1), "Checking 'convergence_checks':\n\tThe 'max_error' must be equal or higher than 0.")
  expect_error(set_convergence_checks(max_SD_error=  -.1), "Checking 'convergence_checks':\n\tThe 'max_SD_error' must be equal or higher than 0.")
  expect_error(set_convergence_checks(max_SD_error=  1.1), "Checking 'convergence_checks':\n\tThe 'max_SD_error' must be equal or lower than 1.")


  expect_equal(set_convergence_checks(), list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(max_Rhat = 1.01),  list(
    max_Rhat      = 1.01,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(min_ESS = 200),  list(
    max_Rhat      = 1.05,
    min_ESS       = 200,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(max_error = 0.01),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = 0.01,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(max_SD_error = 0.01),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = 0.01,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(remove_failed = TRUE),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = TRUE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(balance_probability = FALSE),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = FALSE
  ))

})
