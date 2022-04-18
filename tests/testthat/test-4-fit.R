context("(4) Fitting and updating functions")
skip_on_cran()
skip_on_covr()

# test objects
saved_files <- paste0("fit_", 1:6, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}

# functions simplifying the comparison
remove_time  <- function(fit){
  for(m in 1:length(fit$models)){
    if(is.null(fit$models[[m]]$fit))next
    fit$models[[m]]$fit$timetaken       <- NULL
    fit$models[[m]]$fit$runjags.version <- NULL
  }
  return(fit)
}
clean_all    <- function(fit, only_samples = TRUE){
  if(only_samples){
    fit$data     <- NULL
    fit$add_info <- NULL
    fit$control  <- NULL
    fit$models   <- NULL
  }
  return(fit)
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


test_that("Default model (RoBMA-PSMA) works", {

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
  expect_equal(as.data.frame(fit1$models[[5]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.01)
  expect_equal(as.data.frame(fit1$models[[5]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.01)

  # compare to the previous version object for consistency
  fit1 <- remove_time(fit1)
  expect_equal(clean_all(saved_fits[[1]]), clean_all(fit1))


  ### with predictors, no test
  fit2 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3,
    test_predictors = FALSE,
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 2
  ))

  # compare parameter estimates to flexsurv
  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~  x_cont + x_bin + x_fac3, data = df, dist = "exp")
  expect_equal(as.data.frame(fit2$models[[1]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[1]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.020
  )

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "weibull")
  expect_equal(as.data.frame(fit2$models[[2]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(as.data.frame(fit2$models[[2]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[2]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.020
  )

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "lognorm")
  expect_equal(as.data.frame(fit2$models[[3]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(as.data.frame(fit2$models[[3]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[3]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.020
  )

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "llogis")
  expect_equal(as.data.frame(fit2$models[[4]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(as.data.frame(fit2$models[[4]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[4]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.020
  )

  fitf <- flexsurv::flexsurvreg(Surv(time = time, event = event) ~ x_cont + x_bin + x_fac3, data = df, dist = "gamma")
  expect_equal(as.data.frame(fit2$models[[5]]$fit_summary)["intercept", "Mean"], extract_flexsurv(fitf)[["intercept"]][["mean"]], tolerance = 0.025)
  expect_equal(as.data.frame(fit2$models[[5]]$fit_summary)["aux",       "Mean"], extract_flexsurv(fitf)[["aux"]][["mean"]],       tolerance = 0.025)
  expect_equal(
    as.data.frame(fit2$models[[5]]$fit_summary)[c("x_cont", "x_bin[1]", "x_fac3[B]", "x_fac3[C]"), "Mean"],
    unname(sapply(extract_flexsurv(fitf)$parameters, function(p) p[["mean"]])),
    tolerance = 0.020
  )

  # compare to the previous version object for consistency
  fit2 <- remove_time(fit2)
  expect_equal(clean_all(saved_fits[[2]]), clean_all(fit2))


  ### test continuous predictor
  fit3 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_cont,
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 3
  ))

  fit3 <- remove_time(fit3)
  expect_equal(clean_all(saved_fits[[3]]), clean_all(fit3))


  ### test treatment contrast
  fit4 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_bin,
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 4
  ))

  fit4 <- remove_time(fit4)
  expect_equal(clean_all(saved_fits[[4]]), clean_all(fit4))


  ### test orthonormal contrast
  fit5 <- try_parallel(RoBSA(
    Surv(time = time, event = event) ~ x_fac3,
    priors = list("x_fac3" = prior_factor("mnormal", list(0, 0.25), contrast = "orthonormal")),
    test_predictors = "x_fac3",
    data = df, rescale_data = TRUE, parallel = TRUE, seed = 5
  ))

  fit5 <- remove_time(fit5)
  expect_equal(clean_all(saved_fits[[5]]), clean_all(fit5))


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

  fit6 <- remove_time(fit6)
  expect_equal(clean_all(saved_fits[[6]]), clean_all(fit6))

})

#### creating / updating the test settings ####
if(FALSE){

  saved_fits <- list(fit1, fit2, fit3, fit4, fit5, fit6)

  for(i in 1:length(saved_fits)){
    saved_fits[[i]] <- remove_time(saved_fits[[i]])
  }

  for(i in 1:length(saved_fits)){
    saveRDS(saved_fits[[i]], file = file.path("tests/results/fits/", paste0("fit_",i,".RDS")), compress  = "xz")
  }
}
