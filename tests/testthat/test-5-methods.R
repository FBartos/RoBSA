context("(5) Print and summary functions")
skip_on_cran()

# the summary tables and print functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:6, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}


test_that("Print functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(saved_fits[[i]], print = TRUE, width = 150),
      read.table(file = file.path("../results/print", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }
})

test_that("Summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]]), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # all options
  expect_equal(
    capture_output_lines(summary(saved_fits[[6]], conditional = TRUE, logBF = TRUE, BF01 = TRUE, probs = c(0.10, 0.50, .90)), print = TRUE, width = 150),
    c("Call:"                                                                                       ,
      "RoBSA(formula = Surv(time = time, event = event) ~ x_cont + x_fac3, "                        ,
      "    data = df, priors = list(x_fac3 = list(alt = prior_factor(\"beta\", "                    ,
      "        list(3, 3), contrast = \"treatment\"), null = prior_factor(\"uniform\", "            ,
      "        list(-0.1, 0.1), contrast = \"treatment\"))), test_predictors = \"x_fac3\", "        ,
      "    distributions = c(\"gamma-aft\", \"weibull-aft\", \"lnorm-aft\", "                       ,
      "        \"llogis-aft\"), distributions_weights = c(3, 1, 1, 1), "                            ,
      "    prior_intercept = list(`gamma-aft` = prior(\"normal\", list(1, "                         ,
      "        1)), `weibull-aft` = prior(\"normal\", list(2, 1)), `lnorm-aft` = prior(\"normal\", ",
      "        list(1, 3)), `llogis-aft` = prior(\"normal\", list(1, 4))), "                        ,
      "    parallel = TRUE, seed = 6, rescale_data = TRUE)"                                         ,
      ""                                                                                            ,
      "Robust Bayesian survival analysis"                                                           ,
      "Distributions summary:"                                                                      ,
      "            Models Prior prob. Post. prob. log(Exclusion BF)"                                ,
      "gamma-aft      2/8       0.500       0.914            -2.360"                                ,
      "weibull-aft    2/8       0.167       0.084             0.783"                                ,
      "lnorm-aft      2/8       0.167       0.000            10.029"                                ,
      "llogis-aft     2/8       0.167       0.003             4.370"                                ,
      ""                                                                                            ,
      "Components summary:"                                                                         ,
      "       Models Prior prob. Post. prob. log(Exclusion BF)"                                     ,
      "x_fac3    4/8       0.500       0.189             1.453"                                     ,
      ""                                                                                            ,
      "Model-averaged estimates:"                                                                   ,
      "           Mean Median    0.1   0.5   0.9"                                                   ,
      "x_cont    0.056  0.056 -0.033 0.056 0.146"                                                   ,
      "x_fac3[B] 0.053  0.020 -0.076 0.020 0.275"                                                   ,
      "x_fac3[C] 0.093  0.052 -0.051 0.052 0.376"                                                   ,
      ""                                                                                            ,
      "Conditional estimates:"                                                                      ,
      "           Mean Median    0.1   0.5   0.9"                                                   ,
      "x_cont    0.056  0.056 -0.033 0.056 0.146"                                                   ,
      "x_fac3[B] 0.291  0.284  0.144 0.284 0.444"                                                   ,
      "x_fac3[C] 0.390  0.385  0.228 0.385 0.559"
      )
  )


  # no conditional, yet requested
  expect_equal(
    capture_output_lines(summary(saved_fits[[1]], conditional = TRUE), print = TRUE, width = 150),
    c("Call:"                                                            ,
      "RoBSA(formula = Surv(time = time, event = event) ~ 1, data = df, ",
      "    parallel = TRUE, seed = 1, rescale_data = TRUE)"              ,
      ""                                                                 ,
      "Robust Bayesian survival analysis"                                ,
      "Distributions summary:"                                           ,
      "            Models Prior prob. Post. prob. Inclusion BF"          ,
      "exp-aft        1/5       0.200       0.852       23.061"          ,
      "weibull-aft    1/5       0.200       0.062        0.266"          ,
      "lnorm-aft      1/5       0.200       0.000        0.000"          ,
      "llogis-aft     1/5       0.200       0.001        0.004"          ,
      "gamma-aft      1/5       0.200       0.085        0.370"          ,
      ""                                                                 ,
      "Conditional estimates:"                                           ,
      "[1] Mean   Median 0.025  0.975 "                                  ,
      "<0 rows> (or 0-length row.names)"                                 ,
      ""                                                                 ,
      "Conditional estimates:"                                           ,
      "[1] Mean   Median 0.025  0.975 "                                  ,
      "<0 rows> (or 0-length row.names)"    )
  )


  # distributional parameters & exponential
  expect_equal(
    capture_output_lines(summary(saved_fits[[2]], parameters = TRUE, exp = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                           ,
      "RoBSA(formula = Surv(time = time, event = event) ~ x_cont + x_bin + "            ,
      "    x_fac3, data = df, test_predictors = FALSE, parallel = TRUE, "               ,
      "    seed = 2, rescale_data = TRUE)"                                              ,
      ""                                                                                ,
      "Robust Bayesian survival analysis"                                               ,
      "Distributions summary:"                                                          ,
      "            Models Prior prob. Post. prob. Inclusion BF"                         ,
      "exp-aft        1/5       0.200       0.842       21.310"                         ,
      "weibull-aft    1/5       0.200       0.067        0.286"                         ,
      "lnorm-aft      1/5       0.200       0.000        0.000"                         ,
      "llogis-aft     1/5       0.200       0.004        0.018"                         ,
      "gamma-aft      1/5       0.200       0.087        0.380"                         ,
      ""                                                                                ,
      "Model-averaged estimates:"                                                       ,
      "            Mean Median  0.025 0.975 exp(Mean) exp(Median) exp(0.025) exp(0.975)",
      "x_cont     0.058  0.059 -0.076 0.194     1.060       1.061      0.926      1.214",
      "x_bin[1]  -0.171 -0.173 -0.446 0.098     0.843       0.841      0.640      1.103",
      "x_fac3[B]  0.099  0.099 -0.237 0.431     1.104       1.104      0.789      1.539",
      "x_fac3[C]  0.253  0.250 -0.085 0.592     1.288       1.284      0.919      1.808",
      ""                                                                                ,
      "Distribution estimates (intercept):"                                             ,
      "             Mean Median  0.025 0.975"                                           ,
      "exp-aft     0.521  0.519  0.265 0.769"                                           ,
      "weibull-aft 0.517  0.515  0.252 0.806"                                           ,
      "lnorm-aft   0.027  0.025 -0.303 0.362"                                           ,
      "llogis-aft  0.013  0.009 -0.294 0.345"                                           ,
      "gamma-aft   0.566  0.560  0.246 0.918"                                           ,
      ""                                                                                ,
      "Distribution estimates (auxiliary):"                                             ,
      "             Mean Median 0.025 0.975"                                            ,
      "weibull-aft 0.968  0.967 0.863 1.074"                                            ,
      "lnorm-aft   1.474  1.471 1.336 1.627"                                            ,
      "llogis-aft  1.244  1.243 1.104 1.391"                                            ,
      "gamma-aft   0.969  0.968 0.828 1.120"
      )
  )

  })

test_that("Models summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "models"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary.models", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # test short names
  expect_equal(
    capture_output_lines(summary(saved_fits[[6]], type = "models", short_name = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                              ,
      "RoBSA(formula = Surv(time = time, event = event) ~ x_cont + x_fac3, "                                               ,
      "    data = df, priors = list(x_fac3 = list(alt = prior_factor(\"beta\", "                                           ,
      "        list(3, 3), contrast = \"treatment\"), null = prior_factor(\"uniform\", "                                   ,
      "        list(-0.1, 0.1), contrast = \"treatment\"))), test_predictors = \"x_fac3\", "                               ,
      "    distributions = c(\"gamma-aft\", \"weibull-aft\", \"lnorm-aft\", "                                              ,
      "        \"llogis-aft\"), distributions_weights = c(3, 1, 1, 1), "                                                   ,
      "    prior_intercept = list(`gamma-aft` = prior(\"normal\", list(1, "                                                ,
      "        1)), `weibull-aft` = prior(\"normal\", list(2, 1)), `lnorm-aft` = prior(\"normal\", "                       ,
      "        list(1, 3)), `llogis-aft` = prior(\"normal\", list(1, 4))), "                                               ,
      "    parallel = TRUE, seed = 6, rescale_data = TRUE)"                                                                ,
      ""                                                                                                                   ,
      "Robust Bayesian survival analysis"                                                                                  ,
      "Models overview:"                                                                                                   ,
      " Model Distribution Prior x_cont           Prior x_fac3           Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1    gamma-aft      N(0, 1) treatment contrast: U(-0.1, 0.1)       0.250      -320.81       0.739        8.511",
      "     2  weibull-aft      N(0, 1) treatment contrast: U(-0.1, 0.1)       0.083      -322.07       0.070        0.828",
      "     3    lnorm-aft      N(0, 1) treatment contrast: U(-0.1, 0.1)       0.083      -331.58       0.000        0.000",
      "     4   llogis-aft      N(0, 1) treatment contrast: U(-0.1, 0.1)       0.083      -326.19       0.001        0.013",
      "     5    gamma-aft      N(0, 1)      treatment contrast: B(3, 3)       0.250      -322.26       0.174        0.633",
      "     6  weibull-aft      N(0, 1)      treatment contrast: B(3, 3)       0.083      -323.70       0.014        0.154",
      "     7    lnorm-aft      N(0, 1)      treatment contrast: B(3, 3)       0.083      -331.94       0.000        0.000",
      "     8   llogis-aft      N(0, 1)      treatment contrast: B(3, 3)       0.083      -325.99       0.001        0.015"
    ))

  # test no spikes
  expect_equal(
    capture_output_lines(summary(saved_fits[[5]], type = "models", remove_spike_0 = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                       ,
      "RoBSA(formula = Surv(time = time, event = event) ~ x_fac3, data = df, "                                      ,
      "    priors = list(x_fac3 = prior_factor(\"mnormal\", list(0, 0.25), "                                        ,
      "        contrast = \"orthonormal\")), test_predictors = \"x_fac3\", "                                        ,
      "    parallel = TRUE, seed = 5, rescale_data = TRUE)"                                                         ,
      ""                                                                                                            ,
      "Robust Bayesian survival analysis"                                                                           ,
      "Models overview:"                                                                                            ,
      " Model Distribution              Prior x_fac3              Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1      exp-aft                                              0.100      -317.64       0.601       13.580",
      "     2  weibull-aft                                              0.100      -320.26       0.044        0.412",
      "     3    lnorm-aft                                              0.100      -329.97       0.000        0.000",
      "     4   llogis-aft                                              0.100      -324.44       0.001        0.006",
      "     5    gamma-aft                                              0.100      -319.95       0.060        0.571",
      "     6      exp-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -318.51       0.250        2.996",
      "     7  weibull-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -321.10       0.019        0.172",
      "     8    lnorm-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -330.17       0.000        0.000",
      "     9   llogis-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -324.24       0.001        0.007",
      "    10    gamma-aft orthonormal contrast: mNormal(0, 0.25)       0.100      -320.81       0.025        0.232"
    ))
})

test_that("Diagnostics summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "diagnostics"), print = TRUE, width = 200),
      read.table(file = file.path("../results/summary.diagnostics", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }
})

test_that("Individual summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "individual"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary.individual", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

})


#### creating / updating the test settings ####
if(FALSE){

  saved_files <- paste0("fit_", 1:6, ".RDS")
  saved_fits  <- list()
  for(i in seq_along(saved_files)){
    saved_fits[[i]] <- readRDS(file = file.path("tests/results/fits", saved_files[i]))
  }

  # generate print files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(saved_fits[[i]], print = TRUE, width = 150), file = file.path("tests/results/print", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]]), print = TRUE, width = 150), file = file.path("tests/results/summary", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.models files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "models"), print = TRUE, width = 150), file = file.path("tests/results/summary.models", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.diagnostics files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "diagnostics"), print = TRUE, width = 200), file = file.path("tests/results/summary.diagnostics", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.individual files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "individual"), print = TRUE, width = 150), file = file.path("tests/results/summary.individual", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

}
