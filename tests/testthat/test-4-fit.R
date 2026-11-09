context("(4) Fitting and updating functions")
skip_on_cran()


# test objects
saved_files <- paste0("fit_", 1:6, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}
df <- readRDS(file = file.path("../results/fits", "df.RDS"))

# functions simplifying the comparison
remove_time  <- function(fit){
  for(m in 1:length(fit$models)){
    if(is.null(fit$models[[m]]$fit))next
    fit$models[[m]]$fit$timetaken       <- NULL
    fit$models[[m]]$fit$runjags.version <- NULL
  }
  return(fit)
}
# A refit matches a saved fit only up to Monte Carlo error: compilers, math libraries
# and the state of the R process can change the draws. Each bound is about twice the
# largest difference from the saved fits across 24 other seeds of these six fits.
saved_fit_bounds <- c(log_marglik = 0.02, log_BF = 0.015, post_probs = 0.002, model_means = 0.05, mean = 0.025)
saved_fit_values <- function(fit){
  fit_summary <- summary(fit, conditional = TRUE)
  inference   <- c(fit$RoBSA$inference, fit$RoBSA$inference_distributions)
  values      <- list(
    log_marglik = vapply(fit$models, function(model) model$marglik$logml, numeric(1)),
    log_BF      = vapply(inference, attr, numeric(1), "log_BF"),
    post_probs  = inference[[1]]$post_probs,
    model_means = unlist(lapply(fit$models, function(model) as.data.frame(model$fit_summary)[["Mean"]])),
    mean        = c(as.data.frame(fit_summary$estimates)$Mean, as.data.frame(fit_summary$estimates_conditional)$Mean)
  )
  return(values[lengths(values) > 0])
}
expect_saved_fit <- function(saved_fit, fit){
  saved_values <- saved_fit_values(saved_fit)
  fit_values   <- saved_fit_values(fit)
  expect_identical(lengths(fit_values), lengths(saved_values))
  for(value in names(saved_values)){
    expect_lt(max(abs(fit_values[[value]] - saved_values[[value]])), saved_fit_bounds[[value]], label = value)
  }
}
try_parallel <- function(x, rep = 3){
  temp_fit <- NULL
  i        <- 0
  while(is.null(temp_fit) & i < rep){
    temp_fit <- tryCatch(eval(x), error = function(e) NULL)
    i        <- i + 1
  }
  return(temp_fit)
}


test_that("Fitting function works", {

  ### intercept only
  fit1 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ 1,
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 1
  ))

  # compare parameter estimates to flexsurv
  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ 1, data = df, dist = "exp")
  expect_equal(as.data.frame(fit1$models[[1]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.01)

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ 1, data = df, dist = "weibull")
  expect_equal(as.data.frame(fit1$models[[2]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.01)
  expect_equal(as.data.frame(fit1$models[[2]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.01)

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ 1, data = df, dist = "lognorm")
  expect_equal(as.data.frame(fit1$models[[3]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.01)
  expect_equal(as.data.frame(fit1$models[[3]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.01)

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ 1, data = df, dist = "llogis")
  expect_equal(as.data.frame(fit1$models[[4]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.01)
  expect_equal(as.data.frame(fit1$models[[4]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.01)

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ 1, data = df, dist = "gamma")
  expect_equal(as.data.frame(fit1$models[[5]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.05)
  expect_equal(as.data.frame(fit1$models[[5]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.05)

  # compare to the previous version object for consistency
  expect_saved_fit(saved_fits[[1]], fit1)


  ### with predictors, no test
  fit2 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3,
    test_predictors = FALSE,
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 2
  ))

  # compare parameter estimates to flexsurv (the posterior means differ from the MLEs by
  # the N(0, 1) prior shrinkage, up to about 0.015, and by Monte Carlo error)
  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~  x_cont + x_bin + x_fac3, data = df, dist = "exp")
  expect_equal(as.data.frame(fit2$models[[1]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[1]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.025
  )

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "weibull")
  expect_equal(as.data.frame(fit2$models[[2]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(as.data.frame(fit2$models[[2]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[2]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.025
  )

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "lognorm")
  expect_equal(as.data.frame(fit2$models[[3]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(as.data.frame(fit2$models[[3]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[3]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.025
  )

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "llogis")
  expect_equal(as.data.frame(fit2$models[[4]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(as.data.frame(fit2$models[[4]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[4]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.025
  )

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "gamma")
  expect_equal(as.data.frame(fit2$models[[5]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(as.data.frame(fit2$models[[5]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[5]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.025
  )

  # compare to the previous version object for consistency
  expect_saved_fit(saved_fits[[2]], fit2)


  ### test continuous predictor
  fit3 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_cont,
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 3
  ))

  expect_saved_fit(saved_fits[[3]], fit3)


  ### test treatment contrast
  fit4 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_bin,
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 4
  ))

  expect_saved_fit(saved_fits[[4]], fit4)


  ### test orthonormal contrast
  fit5 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_fac3,
    priors = list("x_fac3" = prior_factor("mnormal", list(0, 0.25), contrast = "orthonormal")),
    test_predictors = "x_fac3",
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 5
  ))

  expect_saved_fit(saved_fits[[5]], fit5)


  ### test more custom model
  fit6 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_cont + x_fac3,
    priors = list("x_fac3" = list(
      "alt"  = prior_factor("beta",    list(3, 3),      contrast = "treatment"),
      "null" = prior_factor("uniform", list(-0.1, 0.1), contrast = "treatment")
    )),
    distributions         = c("gamma-aft", "weibull-aft", "lnorm-aft", "llogis-aft"),
    distributions_weights = c(3, 1, 1, 1),
    prior_intercept = list(
      "gamma-aft"   = prior("normal", list(1, 1)),
      "weibull-aft" = prior("normal", list(2, 1)),
      "lnorm-aft"   = prior("normal", list(1, 3)),
      "llogis-aft"  = prior("normal", list(1, 4))
    ),
    test_predictors = "x_fac3",
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 6
  ))

  expect_saved_fit(saved_fits[[6]], fit6)

})

test_that("Updating function works", {

  fit7 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ 1, distributions = c("exp-aft"),
    data = df, parallel = TRUE, seed = 7
  ))

  # add distribution
  fit7 <- suppressWarnings(update(fit7, distribution = "lnorm-aft"))
  # add a predictor (TODO: specify tests for the added predictor)
  fit7 <- suppressWarnings(update(fit7, formula = Surv(time = time, event = event) ~ x_cont, distribution = "gamma-aft"))
  # change model weights
  fit7 <- suppressWarnings(update(fit7,  model_weights = c(1, 2, 3)))


})

test_that("Parallel fitting works with fewer cores than chains", {

  max_cores <- RoBSA.get_option("max_cores")
  on.exit(RoBSA.options(max_cores = max_cores))
  RoBSA.options(max_cores = 2)

  fit <- RoBSA(
    Surv(time = time, event = event) ~ 1, data = df, distributions = "exp-aft",
    sample = 500, burnin = 200, adapt = 100, parallel = TRUE, seed = 1
  )

  expect_s3_class(fit, "RoBSA")
})

#### creating / updating the test settings ####
if(FALSE){

  # create mock data
  set.seed(68)
  df <- data.frame(
    time   = rgamma(360, 1, 1),
    event  = rbinom(360, 1, 0.5),
    x_cont = rnorm(360),
    x_bin  = as.factor(rbinom(360, 1, .5)),
    x_fac3 = factor(rep(c("A", "B", "C"), 120), levels = c("A", "B", "C")),
    x_fac4 = factor(rep(c("A", "B", "C", "D"), 90), levels = c("A", "B", "C", "D"))
  )
  saveRDS(df, file = file.path("tests/results/fits/", "df.RDS"), compress  = "xz")


  saved_fits <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7)

  for(i in 1:length(saved_fits)){
    saved_fits[[i]] <- remove_time(saved_fits[[i]])
  }

  for(i in 1:length(saved_fits)){
    saveRDS(saved_fits[[i]], file = file.path("tests/results/fits/", paste0("fit_",i,".RDS")), compress  = "xz")
  }
}
