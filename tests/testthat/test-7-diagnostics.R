context("(7) Diagnostics plots")
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

test_that("Diagnostics plots work", {

  # treatment contrast
  expect_doppelganger(paste0("diagnostics_autocorrelation_1"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(4, 2))
    suppressMessages(diagnostics_autocorrelation(saved_fits[[6]], "x_fac3"))
  })
  expect_doppelganger(paste0("diagnostics_autocorrelation_2"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(2, 2))
    suppressMessages(diagnostics_autocorrelation(saved_fits[[6]], "x_fac3", lags = 20, show_models = 7:8))
  })
  expect_doppelganger(paste0("diagnostics_trace_1"),   function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(4, 2))
    suppressMessages(diagnostics_trace(saved_fits[[6]], "x_fac3"))
  })
  expect_doppelganger(paste0("diagnostics_density_1"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2))
    diagnostics_density(saved_fits[[6]], "x_fac3", show_models = 8)
  })
  expect_doppelganger(paste0("diagnostics_density_2"), suppressMessages(diagnostics_density(saved_fits[[6]], "x_fac3", plot_type = "ggplot", show_models = 8)[[8]][[1]]))

  # orthonormal contrast
  expect_doppelganger(paste0("diagnostics_autocorrelation_3"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(5, 3))
    suppressMessages(diagnostics_autocorrelation(saved_fits[[5]], "x_fac3"))
  })
  expect_doppelganger(paste0("diagnostics_autocorrelation_4"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(2, 3))
    suppressMessages(diagnostics_autocorrelation(saved_fits[[5]], "x_fac3", lags = 20, show_models = 7:8))
  })
  expect_doppelganger(paste0("diagnostics_trace_2"),   function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 3))
    suppressMessages(diagnostics_trace(saved_fits[[5]], "x_fac3", show_models = 7:9))
  })
  expect_doppelganger(paste0("diagnostics_density_3"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(1, 3))
    diagnostics_density(saved_fits[[5]], "x_fac3", show_models = 8)
  })
  expect_doppelganger(paste0("diagnostics_density_4"), suppressMessages(diagnostics_density(saved_fits[[5]], "x_fac3", plot_type = "ggplot", show_models = 8)[[8]][[2]]))

  # continuous predictor
  expect_doppelganger(paste0("diagnostics_autocorrelation_5"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 2))
    suppressMessages(diagnostics_autocorrelation(saved_fits[[3]], "x_cont"))
  })
  expect_doppelganger(paste0("diagnostics_autocorrelation_6"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2))
    suppressMessages(diagnostics_autocorrelation(saved_fits[[3]], "x_cont", lags = 20, show_models = 7:8))
  })
  expect_doppelganger(paste0("diagnostics_trace_3"),   function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(3, 2))
    suppressMessages(diagnostics_trace(saved_fits[[3]], "x_cont"))
  })
  expect_doppelganger(paste0("diagnostics_density_5"), function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mfrow = c(1, 1))
    diagnostics_density(saved_fits[[3]], "x_cont", show_models = 8, xlim = c(-1, 1))
  })
  expect_doppelganger(paste0("diagnostics_density_6"), suppressMessages(diagnostics_density(saved_fits[[3]], "x_cont", plot_type = "ggplot", show_models = 8)[[8]]))


  # test errors
  expect_error(diagnostics_density(saved_fits[[5]], "x_fac3o"), "The passed parameter does not correspond to any of the specified predictors: 'x_fac3'", fixed = TRUE)
})
