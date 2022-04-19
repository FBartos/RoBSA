context("(6) Plot functions")
skip_on_cran()

# the plotting functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:6, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}
df <- readRDS(file = file.path("../results/fits", "df.RDS"))

set.seed(1)

test_that("Parameter plots work", {

  # treatment contrast
  expect_doppelganger(paste0("plot_est_6_2.1"), function()plot(saved_fits[[6]], "x_fac3"))
  expect_doppelganger(paste0("plot_est_6_2.2"), function()plot(saved_fits[[6]], "x_fac3", prior = TRUE))
  expect_doppelganger(paste0("plot_est_6_2.3"), function()plot(saved_fits[[6]], "x_fac3", conditional = TRUE))
  expect_doppelganger(paste0("plot_est_6_2.4"), function()plot(saved_fits[[6]], "x_fac3", conditional = TRUE, prior = TRUE))

  expect_doppelganger(paste0("ggplot_est_6_2.1"), plot(saved_fits[[6]], plot_type = "ggplot", "x_fac3"))
  expect_doppelganger(paste0("ggplot_est_6_2.2"), plot(saved_fits[[6]], plot_type = "ggplot", "x_fac3", prior = TRUE))
  expect_doppelganger(paste0("ggplot_est_6_2.3"), plot(saved_fits[[6]], plot_type = "ggplot", "x_fac3", conditional = TRUE))
  expect_doppelganger(paste0("ggplot_est_6_2.4"), plot(saved_fits[[6]], plot_type = "ggplot", "x_fac3", conditional = TRUE, prior = TRUE))

  # orthonormal contrast
  expect_doppelganger(paste0("plot_est_5_1.1"), function()plot(saved_fits[[5]], "x_fac3"))
  expect_doppelganger(paste0("plot_est_5_1.2"), function()plot(saved_fits[[5]], "x_fac3", prior = TRUE))
  expect_doppelganger(paste0("plot_est_5_1.3"), function()plot(saved_fits[[5]], "x_fac3", conditional = TRUE))
  expect_doppelganger(paste0("plot_est_5_1.4"), function()plot(saved_fits[[5]], "x_fac3", conditional = TRUE, prior = TRUE))

  expect_doppelganger(paste0("ggplot_est_5_1.1"), plot(saved_fits[[5]], plot_type = "ggplot", "x_fac3"))
  expect_doppelganger(paste0("ggplot_est_5_1.2"), plot(saved_fits[[5]], plot_type = "ggplot", "x_fac3", prior = TRUE))
  expect_doppelganger(paste0("ggplot_est_5_1.3"), plot(saved_fits[[5]], plot_type = "ggplot", "x_fac3", conditional = TRUE))
  expect_doppelganger(paste0("ggplot_est_5_1.4"), plot(saved_fits[[5]], plot_type = "ggplot", "x_fac3", conditional = TRUE, prior = TRUE))

  # continuous predictor
  expect_doppelganger(paste0("plot_est_3_1.1"), function()plot(saved_fits[[3]], "x_cont"))
  expect_doppelganger(paste0("plot_est_3_1.2"), function()plot(saved_fits[[3]], "x_cont", prior = TRUE))
  expect_doppelganger(paste0("plot_est_3_1.3"), function()plot(saved_fits[[3]], "x_cont", conditional = TRUE))
  expect_doppelganger(paste0("plot_est_3_1.4"), function()plot(saved_fits[[3]], "x_cont", conditional = TRUE, prior = TRUE))

  expect_doppelganger(paste0("ggplot_est_3_1.1"), plot(saved_fits[[3]], plot_type = "ggplot", "x_cont"))
  expect_doppelganger(paste0("ggplot_est_3_1.2"), plot(saved_fits[[3]], plot_type = "ggplot", "x_cont", prior = TRUE))
  expect_doppelganger(paste0("ggplot_est_3_1.3"), plot(saved_fits[[3]], plot_type = "ggplot", "x_cont", conditional = TRUE))
  expect_doppelganger(paste0("ggplot_est_3_1.4"), plot(saved_fits[[3]], plot_type = "ggplot", "x_cont", conditional = TRUE, prior = TRUE))

  # additional settings
  expect_doppelganger(paste0("plot_est_3_1.5"), function()plot(saved_fits[[3]], "x_cont", prior = TRUE, dots_prior = list(col = "blue", lty = 2), col = "red", lty = 2, xlim = c(0, 1), main = "Title"))
  expect_doppelganger(paste0("ggplot_est_5_1.6"), function()plot(saved_fits[[5]], "x_fac3", prior = TRUE, dots_prior = list(col = "grey", lty = 2), col = c("blue", "red", "orange"),  xlim = c(0, 1), main = "Title"))
  expect_doppelganger(paste0("ggplot_est_5_1.7"), function()plot(saved_fits[[6]], "x_fac3", prior = TRUE, dots_prior = list(col = "grey", lty = 2), col = c("blue", "green"),  xlim = c(0, 1), main = "Title"))

  # plot intercepts & auxiliary parameters
  expect_doppelganger(paste0("ggplot_est_5_1"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 2))
    plot(saved_fits[[1]], "intercept", prior = TRUE)
  })
  expect_doppelganger(paste0("ggplot_est_5_2"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 2))
    plot(saved_fits[[1]], "aux", prior = TRUE)
  })

  # test errors
  expect_error(plot(saved_fits[[5]], "x_fac3o"), "The passed parameter does not correspond to any of the specified predictors: 'x_fac3'", fixed = TRUE)
})


test_that("Individual model plots work", {

  # treatment contrast
  expect_doppelganger(paste0("plot_models_est_6_2.1"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 1))
    plot_models(saved_fits[[6]], parameter = "x_fac3")
  })
  expect_doppelganger(paste0("plot_models_est_6_2.2"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 1))
    plot_models(saved_fits[[6]], parameter = "x_fac3", conditional = TRUE)
  })

  # orthonormal contrast
  expect_doppelganger(paste0("plot_models_est_5_1.1"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 1))
    plot_models(saved_fits[[5]], parameter = "x_fac3")
  })
  expect_doppelganger(paste0("plot_models_est_5_1.2"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 1))
    plot_models(saved_fits[[5]], parameter = "x_fac3", conditional = TRUE)
  })

  # continuous predictor
  expect_doppelganger(paste0("plot_models_3_1.1"), plot_models(saved_fits[[3]], parameter = "x_cont", plot_type = "ggplot"))
  expect_doppelganger(paste0("plot_models_3_1.2"), plot_models(saved_fits[[3]], parameter = "x_cont", plot_type = "ggplot", conditional = TRUE))

  # test errors
  expect_error(plot(saved_fits[[3]], "x_fac3"), "The passed parameter does not correspond to any of the specified predictors: 'x_cont'", fixed = TRUE)
})


test_that("Survival, hazard, and density plots work", {

  # compare with flexsurv
  fit2_exp     <- update(saved_fits[[2]], model_weights = c(1, 0, 0, 0, 0))
  fit2_weibull <- update(saved_fits[[2]], model_weights = c(0, 1, 0, 0, 0))
  fit2_lnorm   <- update(saved_fits[[2]], model_weights = c(0, 0, 1, 0, 0))
  fit2_llogis  <- update(saved_fits[[2]], model_weights = c(0, 0, 0, 1, 0))
  fit2_gamma   <- update(saved_fits[[2]], model_weights = c(0, 0, 0, 0, 1))

  fitf_exp     <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "exp")
  fitf_weibull <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "weibull")
  fitf_lnorm   <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "lognorm")
  fitf_llogis  <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "llogis")
  fitf_gamma   <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "gamma")

  expect_doppelganger(paste0("plot_sur_2.3.1"), function(){
    plot_survival(fit2_exp, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_exp, col = "red", ci = FALSE, type = "survival", lwd = 1)
  })

  expect_doppelganger(paste0("plot_sur_2.3.2"), function(){
    plot_survival(fit2_weibull, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_weibull, col = "red", ci = FALSE, type = "survival", lwd = 1)
  })

  expect_doppelganger(paste0("plot_sur_2.3.3"), function(){
    plot_survival(fit2_lnorm, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_lnorm, col = "red", ci = FALSE, type = "survival", lwd = 1)
  })

  expect_doppelganger(paste0("plot_sur_2.3.4"), function(){
    plot_survival(fit2_llogis, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_llogis, col = "red", ci = FALSE, type = "survival", lwd = 1)
  })

  expect_doppelganger(paste0("plot_sur_2.3.5"), function(){
    plot_survival(fit2_gamma, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_gamma, col = "red", ci = FALSE, type = "survival", lwd = 1)
  })

  expect_doppelganger(paste0("plot_haz_2.3.1"), function(){
    plot_hazard(fit2_exp, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_exp, col = "red", ci = FALSE, type = "hazard", lwd = 1)
  })

  expect_doppelganger(paste0("plot_haz_2.3.2"), function(){
    plot_hazard(fit2_weibull, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_weibull, col = "red", ci = FALSE, type = "hazard", lwd = 1)
  })

  expect_doppelganger(paste0("plot_haz_2.3.3"), function(){
    plot_hazard(fit2_lnorm, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_lnorm, col = "red", ci = FALSE, type = "hazard", lwd = 1)
  })

  expect_doppelganger(paste0("plot_haz_2.3.4"), function(){
    plot_hazard(fit2_llogis, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_llogis, col = "red", ci = FALSE, type = "hazard", lwd = 1)
  })

  expect_doppelganger(paste0("plot_haz_2.3.5"), function(){
    plot_hazard(fit2_gamma, new_data = data.frame(x_fac3 = "A", x_bin = 0, x_cont = mean(df$time)), plot_type = "base")
    lines(fitf_gamma, col = "red", ci = FALSE, type = "hazard", lwd = 1)
  })


  # treatment contrast
  expect_doppelganger(paste0("ggplot_sur_6_1.1"), plot_survival(saved_fits[[6]], plot_type = "ggplot", predictor = "x_fac3"))
  expect_doppelganger(paste0("ggplot_haz_6_1.1"), plot_hazard(saved_fits[[6]],   plot_type = "ggplot",   predictor = "x_fac3"))
  expect_doppelganger(paste0("ggplot_den_6_1.1"), plot_density(saved_fits[[6]],  plot_type = "ggplot",  predictor = "x_fac3"))

  # orthonormal contrast
  expect_doppelganger(paste0("plot_sur_5_1.1"), function()plot_survival(saved_fits[[5]], predictor = "x_fac3"))
  expect_doppelganger(paste0("plot_haz_5_1.1"), function()plot_hazard(saved_fits[[5]],   predictor = "x_fac3"))
  expect_doppelganger(paste0("plot_den_5_1.1"), function()plot_density(saved_fits[[5]],  predictor = "x_fac3"))

  # continuous predictor
  expect_doppelganger(paste0("ggplot_sur_3_1.1"), plot_survival(saved_fits[[3]], plot_type = "ggplot", predictor = "x_cont"))
  expect_doppelganger(paste0("ggplot_haz_3_1.2"), plot_hazard(saved_fits[[3]],   plot_type = "ggplot", predictor = "x_cont"))
  expect_doppelganger(paste0("ggplot_den_3_1.3"), plot_density(saved_fits[[3]],  plot_type = "ggplot", predictor = "x_cont", conditional = TRUE))

})


