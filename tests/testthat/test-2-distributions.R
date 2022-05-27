context("(2) Distribution functions")
skip_on_cran()

### weighted normal distributions ----
test_that("Input checks work", {

  expect_equal(RoBSA.get_option("distributions"), c("exp-aft","weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"))

  # check that only times > 0 are supported
  for(distribution in RoBSA.get_option("distributions")){
    for(type in c("log_density", "log_hazard", "log_survival", "density", "hazard", "survival")){
      fun <- eval(parse(text = paste0(gsub("-", "_", distribution), "_", type)))
      if(distribution == "exp-aft"){
        expect_error(fun(t = -rgamma(1, 1, 1), eta = rnorm(1)), "The 't' must be equal or higher than 0.")
      }else{
        expect_error(do.call(fun, list(t = -rgamma(1, 1, 1), eta = rnorm(1), 1)), "The 't' must be equal or higher than 0.")
      }
    }
  }

  # check warnings for undefined moments
  expect_warning(llogis_aft_mean(0.5, 1),              "llogis-aft: mean is undefinied for shape <= 1")
  expect_warning(llogis_aft_mean(0.5, c(0.5, 0.8)),    "llogis-aft: mean is undefinied for shape <= 1")
  expect_warning(llogis_aft_sd(0.5, 2),                "llogis-aft: sd is undefinied for shape <= 2")
  expect_warning(llogis_aft_sd(0.5, c(0.5, 0.8, 1.2)), "llogis-aft: sd is undefinied for shape <= 2")
})

test_that("check consistency between JAGS and R versions", {

  # re-load the module
  RoBSA:::.load_RoBSA_module()

  set.seed(42)
  df <- data.frame(
    time = rgamma(10, 1, 1)
  )

  for(distribution in RoBSA.get_option("distributions")){

    set.seed(1)
    if(distribution == "exp-aft"){

      model_syntax <- paste0(
      "model{\n",
      "  eta  ~ dnorm(0, 1)\n",
      "  time ~ dgamma(1, 1)\n",
      paste0("  log_den = ", gsub("-", "_", distribution), "_event_lpdf(time, eta)\n"),
      paste0("  log_sur = ", gsub("-", "_", distribution), "_cens_r_lpdf(time, eta)\n"),
      "}"
      )

      model <- rjags::jags.model(file = textConnection(model_syntax), quiet = TRUE)
      fit   <- rjags::coda.samples(model = model, variable.names = c("eta", "time", "log_den", "log_sur"), n.iter = 100, quiet = TRUE, progress.bar = "none")

      expect_equal(as.vector(fit[[1]][,"log_den"]), do.call(
        eval(parse(text = paste0(gsub("-", "_", distribution), "_log_density"))),
        list(
          t   = as.vector(fit[[1]][,"time"]),
          eta = as.vector(fit[[1]][,"eta"])
      )))
      expect_equal(as.vector(fit[[1]][,"log_sur"]), do.call(
        eval(parse(text = paste0(gsub("-", "_", distribution), "_log_survival"))),
        list(
          t   = as.vector(fit[[1]][,"time"]),
          eta = as.vector(fit[[1]][,"eta"])
      )))

    }else{

      model_syntax <- paste0(
        "model{\n",
        "  eta  ~ dnorm(0, 1)\n",
        "  aux  ~ dgamma(1, 1)\n",
        "  time ~ dgamma(1, 1)\n",
        paste0("  log_den = ", gsub("-", "_", distribution), "_event_lpdf(time, eta, aux)\n"),
        paste0("  log_sur = ", gsub("-", "_", distribution), "_cens_r_lpdf(time, eta, aux)\n"),
        "}"
      )

      model <- rjags::jags.model(file = textConnection(model_syntax), quiet = TRUE)
      fit   <- rjags::coda.samples(model = model, variable.names = c("eta", "aux", "time", "log_den", "log_sur"), n.iter = 100, quiet = TRUE, progress.bar = "none")

      expect_equal(as.vector(fit[[1]][,"log_den"]), do.call(
        eval(parse(text = paste0(gsub("-", "_", distribution), "_log_density"))),
        list(
          t   = as.vector(fit[[1]][,"time"]),
          eta = as.vector(fit[[1]][,"eta"]),
          as.vector(fit[[1]][,"aux"])
      )))
      expect_equal(as.vector(fit[[1]][,"log_sur"]), do.call(
        eval(parse(text = paste0(gsub("-", "_", distribution), "_log_survival"))),
        list(
          t   = as.vector(fit[[1]][,"time"]),
          eta = as.vector(fit[[1]][,"eta"]),
          as.vector(fit[[1]][,"aux"])
      )))
    }

  }
})

test_that("check consistency between R's pdf, rng, quantile, and momemt functions", {

  for(distribution in RoBSA.get_option("distributions")){

    set.seed(1)
    if(distribution == "exp-aft"){

      expect_doppelganger(distribution, function(){

        dist_fname <- gsub("-", "_", distribution)
        samples    <- do.call(
          eval(parse(text = paste0(dist_fname, "_r"))),
          list(
            n   = 10000,
            eta = 1.5
        ))

        hist(samples, breaks = 100, freq = FALSE, main = distribution)
        curve(eval(parse(text = paste0(dist_fname, "_density")))(x, eta = 1.5), from = 0, to = max(samples), col = "blue", lwd = 2, add = TRUE)

        abline(v = mean(samples), lwd = 2)
        abline(v = eval(parse(text = paste0(dist_fname, "_mean")))(eta = 1.5), lwd = 2, col = "blue", lty = 2)

        abline(v = mean(samples) + sd(samples), lwd = 2)
        abline(v = eval(parse(text = paste0(dist_fname, "_mean")))(eta = 1.5) + eval(parse(text = paste0(dist_fname, "_sd")))(eta = 1.5), lwd = 2, col = "blue", lty = 2)

        abline(v = quantile(samples, probs = 0.10), lwd = 2)
        abline(v = eval(parse(text = paste0(dist_fname, "_q")))(p = 0.10, eta = 1.5), lwd = 2, col = "blue", lty = 2)

        expect_equal(
          mean(samples < mean(samples)),
          eval(parse(text = paste0(dist_fname, "_p")))(q = mean(samples), eta = 1.5),
          tolerance = 1e-2
        )
      })

    }else{

      expect_doppelganger(distribution, function(){

        dist_fname <- gsub("-", "_", distribution)
        samples    <- do.call(
          eval(parse(text = paste0(dist_fname, "_r"))),
          list(
            n   = 10000,
            eta = 5.5,
            2.5
          ))

        hist(samples, breaks = 100, freq = FALSE, main = distribution)
        curve(eval(parse(text = paste0(dist_fname, "_density")))(x, eta = 5.5, 2.5), from = 0, to = max(samples), col = "blue", lwd = 2, add = TRUE)

        abline(v = mean(samples), lwd = 2)
        abline(v = eval(parse(text = paste0(dist_fname, "_mean")))(eta = 5.5, 2.5), lwd = 2, col = "blue", lty = 2)

        abline(v = mean(samples) + sd(samples), lwd = 2)
        abline(v = eval(parse(text = paste0(dist_fname, "_mean")))(eta = 5.5, 2.5) + eval(parse(text = paste0(dist_fname, "_sd")))(eta = 5.5, 2.5), lwd = 2, col = "blue", lty = 2)

        abline(v = quantile(samples, probs = 0.10), lwd = 2)
        abline(v = eval(parse(text = paste0(dist_fname, "_q")))(p = 0.10, eta = 5.5, 2.5), lwd = 2, col = "blue", lty = 2)

        expect_equal(
          mean(samples < mean(samples)),
          eval(parse(text = paste0(dist_fname, "_p")))(q = mean(samples), eta = 5.5, 2.5),
          tolerance = 1e-2
        )
      })

    }

  }
})
