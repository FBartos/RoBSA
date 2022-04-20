context("(8) Predict")
skip_on_cran()

# hazard, density, and survival is already tested in the plots
saved_files <- paste0("fit_", 1:6, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}
df <- readRDS(file = file.path("../results/fits", "df.RDS"))

set.seed(1)

test_that("Predict mean survival works", {

  # test against flexsurf
  prediction_models <- suppressWarnings(predict(saved_fits[[2]], predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "mean", averaged = FALSE))

  fitf_exp     <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "exp")
  fitf_weibull <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "weibull")
  fitf_lnorm   <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "lognorm")
  fitf_llogis  <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "llogis")
  fitf_gamma   <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "gamma")

  predf_exp     <- data.frame(predict(fitf_exp,     data.frame(x_fac3 = c("A", "B", "C"), x_bin = as.factor(0), x_cont = mean(df$x_cont)), type = "response"))[,1]
  predf_weibull <- data.frame(predict(fitf_weibull, data.frame(x_fac3 = c("A", "B", "C"), x_bin = as.factor(0), x_cont = mean(df$x_cont)), type = "response"))[,1]
  predf_lnorm   <- data.frame(predict(fitf_lnorm,   data.frame(x_fac3 = c("A", "B", "C"), x_bin = as.factor(0), x_cont = mean(df$x_cont)), type = "response"))[,1]
  predf_llogis  <- data.frame(predict(fitf_llogis,  data.frame(x_fac3 = c("A", "B", "C"), x_bin = as.factor(0), x_cont = mean(df$x_cont)), type = "response"))[,1]
  predf_gamma   <- data.frame(predict(fitf_gamma,   data.frame(x_fac3 = c("A", "B", "C"), x_bin = as.factor(0), x_cont = mean(df$x_cont)), type = "response"))[,1]

  expect_equal(prediction_models[[1]][,"median"] / predf_exp,     rep(1, 3), tolerance = 0.02)
  expect_equal(prediction_models[[2]][,"median"] / predf_weibull, rep(1, 3), tolerance = 0.02)
  expect_equal(prediction_models[[3]][,"median"] / predf_lnorm,   rep(1, 3), tolerance = 0.03)
  expect_equal(prediction_models[[4]][,"median"] / predf_llogis,  rep(1, 3), tolerance = 0.10)
  # the mean of llogistic distribution is super sensitive to the parameter estimates
  expect_equal(prediction_models[[5]][,"median"] / predf_gamma,   rep(1, 3), tolerance = 0.02)

  # check posterior of the predictions for the llogis distribution
  prediction_models <- suppressWarnings(predict(saved_fits[[2]], predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "mean", averaged = FALSE, summarize = FALSE))
  expect_doppelganger("mean_llogis", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 1))

    hist(prediction_models[[4]][1,], breaks = 50)
    abline(v = predf_llogis[1], lwd = 2, col = "blue")

    hist(prediction_models[[4]][2,], breaks = 50)
    abline(v = predf_llogis[2], lwd = 2, col = "blue")

    hist(prediction_models[[4]][3,], breaks = 50)
    abline(v = predf_llogis[3], lwd = 2, col = "blue")
  })


  # test consistency for SD predictions
  fit2_exp     <- update(saved_fits[[2]], model_weights = c(1, 0, 0, 0, 0))
  fit2_weibull <- update(saved_fits[[2]], model_weights = c(0, 1, 0, 0, 0))
  fit2_lnorm   <- update(saved_fits[[2]], model_weights = c(0, 0, 1, 0, 0))
  fit2_llogis  <- update(saved_fits[[2]], model_weights = c(0, 0, 0, 1, 0))
  fit2_gamma   <- update(saved_fits[[2]], model_weights = c(0, 0, 0, 0, 1))

  prediction_models <- suppressWarnings(predict(saved_fits[[2]], predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "sd", averaged = FALSE))
  expect_equal(as.matrix(prediction_models[[1]]), as.matrix(predict(fit2_exp,     predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "sd", averaged = TRUE)), tolerance = 0.01)
  expect_equal(as.matrix(prediction_models[[2]]), as.matrix(predict(fit2_weibull, predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "sd", averaged = TRUE)), tolerance = 0.01)
  expect_equal(as.matrix(prediction_models[[3]]), as.matrix(predict(fit2_lnorm,   predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "sd", averaged = TRUE)), tolerance = 0.01)
  expect_equal(as.matrix(prediction_models[[4]]), as.matrix(suppressWarnings(predict(fit2_llogis,  predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "sd", averaged = TRUE))), tolerance = 0.01)
  expect_equal(as.matrix(prediction_models[[5]]), as.matrix(predict(fit2_gamma,   predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "sd", averaged = TRUE)), tolerance = 0.01)

  # check posterior of survival prediction
  prediction_models <- suppressWarnings(predict(saved_fits[[2]],t = 0.5,  predictor = "x_fac3", covariates_data = data.frame(x_bin = 0, x_cont = mean(df$x_cont)), type = "survival", averaged = FALSE, summarize = FALSE))
  predf_llogis  <- data.frame(predict(fitf_llogis, times = 0.5,  data.frame(x_fac3 = c("A", "B", "C"), x_bin = as.factor(0), x_cont = mean(df$x_cont)), type = "survival"))[,2]

  expect_doppelganger("surv_llogis", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 1))

    hist(prediction_models[[4]][[1]][1,], breaks = 50)
    abline(v = predf_llogis[1], lwd = 2, col = "blue")

    hist(prediction_models[[4]][[2]][1,], breaks = 50)
    abline(v = predf_llogis[2], lwd = 2, col = "blue")

    hist(prediction_models[[4]][[3]][1,], breaks = 50)
    abline(v = predf_llogis[3], lwd = 2, col = "blue")
  })


})

