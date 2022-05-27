context("(9) Prior calibration functions")
skip_on_cran()


test_that("Quantile calibration works", {

  median_t   <- 5
  iq_range_t <- 3

  out <- calibrate_quartiles(median_t = median_t, iq_range_t = iq_range_t, prior_sd = 0.5, distributions = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"),
                             verbose = FALSE, search_bounds1 = c(-1e2, 1e2), search_bounds2 = c(0 + 1e-2, 1e2))

  expect_equal(exp_aft_q    (0.50, mean(out$intercept$`exp-aft`)), median_t, tolerance = 1e-3)
  expect_equal(weibull_aft_q(0.50, mean(out$intercept$`weibull-aft`), mean(out$aux$`weibull-aft`)), median_t, tolerance = 1e-3)
  expect_equal(lnorm_aft_q  (0.50, mean(out$intercept$`lnorm-aft`),   mean(out$aux$`lnorm-aft`)),   median_t, tolerance = 1e-3)
  expect_equal(llogis_aft_q (0.50, mean(out$intercept$`llogis-aft`),  mean(out$aux$`llogis-aft`)),  median_t, tolerance = 1e-3)
  expect_equal(gamma_aft_q  (0.50, mean(out$intercept$`gamma-aft`),   mean(out$aux$`gamma-aft`)),   median_t, tolerance = 1e-3)

  expect_equal(
    weibull_aft_q(0.75, mean(out$intercept$`weibull-aft`), mean(out$aux$`weibull-aft`))
    - weibull_aft_q(0.25, mean(out$intercept$`weibull-aft`), mean(out$aux$`weibull-aft`)),
    iq_range_t, tolerance = 1e-3)
  expect_equal(
    lnorm_aft_q(0.75, mean(out$intercept$`lnorm-aft`),   mean(out$aux$`lnorm-aft`))
    - lnorm_aft_q(0.25, mean(out$intercept$`lnorm-aft`),   mean(out$aux$`lnorm-aft`)),
    iq_range_t, tolerance = 1e-3)
  expect_equal(
    llogis_aft_q(0.75, mean(out$intercept$`llogis-aft`),  mean(out$aux$`llogis-aft`))
    - llogis_aft_q(0.25, mean(out$intercept$`llogis-aft`),  mean(out$aux$`llogis-aft`)),
    iq_range_t, tolerance = 1e-3)
  expect_equal(
    gamma_aft_q(0.75, mean(out$intercept$`gamma-aft`),   mean(out$aux$`gamma-aft`))
    -gamma_aft_q(0.25, mean(out$intercept$`gamma-aft`),   mean(out$aux$`gamma-aft`)),
    iq_range_t, tolerance = 1e-3)

})

test_that("Meta-analytic calibration works", {

  set.seed(89)
  dfs <- list(
    data.frame(
      time   = gamma_aft_r(50, 0.5, 1.5),
      status = stats::rbinom(50, 1, 0.5)
    ),
    data.frame(
      time   = gamma_aft_r(150, 0.5, 1.5),
      status = stats::rbinom(150, 1, 0.5)
    ),
    data.frame(
      time   = gamma_aft_r(70, 0.5, 1.5),
      status = stats::rbinom(70, 1, 0.5)
    ),
    data.frame(
      time   = gamma_aft_r(60, 0.5, 1.5),
      status = stats::rbinom(60, 1, 0.5)
    ),
    data.frame(
      time   = gamma_aft_r(50, 0.5, 1.5),
      status = stats::rbinom(50, 1, 0.5)
    )
  )


  out <- calibrate_meta_analytic(dfs)

  expect_doppelganger("calibrate_meta_analytic", function(){
    curve(gamma_aft_density(t = x, eta = 0.5, shape = 1.5), from = 0, to = 15, lwd = 2, las = 1, xlab = "Time", ylab = "Density", bty = "n")
    curve(exp_aft_density(t = x, eta = mean(out$intercept$`exp-aft`)), from = 0, to = 15, add = TRUE)
    curve(weibull_aft_density(t = x, eta = mean(out$intercept$`weibull-aft`), shape = mean(out$aux$`weibull-aft`)), from = 0, to = 15, add = TRUE)
    curve(lnorm_aft_density(t = x, eta = mean(out$intercept$`lnorm-aft`), sd = mean(out$aux$`lnorm-aft`)), from = 0, to = 15, add = TRUE)
    curve(llogis_aft_density(t = x, eta = mean(out$intercept$`llogis-aft`), shape = mean(out$aux$`llogis-aft`)), from = 0, to = 15, add = TRUE)
    curve(gamma_aft_density(t = x, eta = mean(out$intercept$`gamma-aft`), shape = mean(out$aux$`gamma-aft`)), from = 0, to = 15, add = TRUE)
  })

})
