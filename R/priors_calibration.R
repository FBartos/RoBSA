#' @title Calibrate prior distributions based on quartiles
#'
#' @description Calibrates prior distributions for parametric
#' survival analysis based on median survival and interquartile
#' range. Returns a list of prior distribution for the intercepts
#' and auxiliary parameters.
#'
#' @param median_t median survival
#' @param iq_range_t interquartile range of the survival
#' @param prior_sd pre-specified standard deviation of the
#' prior distributions (either a single value that is used
#' for both the intercept and auxiliary parameter or a vector
#' where the first value corresponds to the sd for the prior
#' distribution on the intercept and the second value to the
#' sd for the prior distribution on the auxiliary parameter)
#' @param distributions vector of parametric families for which
#' prior distributions ought to be calibrated
#' @param verbose whether debug information be printed
#' @param search_bounds1 search boundaries for the intercept
#' parameter
#' @param search_bounds2 search boundaries for the auxiliary
#' parameter
#'
#' @export
calibrate_quartiles      <- function(median_t, iq_range_t, prior_sd = 0.5, distributions = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"),
                                     verbose = FALSE, search_bounds1 = c(-1e2, 1e2), search_bounds2 = c(0 + 1e-2, 1e2)){

  parameters <- list()
  priors     <- list()
  priors_int <- list()
  priors_aux <- list()

  # set the same sd for both parameters if unspecified
  if(length(prior_sd) == 1){
    prior_sd <- rep(prior_sd, 2)
  }

  for(distribution in distributions){

    parameters[[distribution]] <- .calibrate_quantiles_fun(median_t, iq_range_t, distribution, verbose, search_bounds1, search_bounds2)
    priors[[distribution]]     <- .make_prior(parameters[[distribution]], prior_sd, distribution)

    priors_int[[distribution]] <- priors[[distribution]][["intercept"]]
    if(.has_aux(distribution)){
      priors_aux[[distribution]] <- priors[[distribution]][["aux"]]
    }
  }

  return(list(
    intercept = priors_int,
    aux       = priors_aux
  ))
}
.calibrate_quantiles_fun <- function(median_t, iq_range_t, distribution, verbose, search_bounds1, search_bounds2){

  if(verbose){
    cat(paste0("Distribution: ", distribution, "\n"))
  }


  if(.has_aux(distribution)){

    parameters <- optim(
      par   = c(0, 1),
      fn    = .solve_quantiles2,
      lower = c(search_bounds1[1], search_bounds2[1]),
      upper = c(search_bounds1[2], search_bounds2[2]),
      median_t     = median_t,
      iq_range_t   = iq_range_t,
      quant_fun    = eval(parse(text = paste0(".", gsub("-", "_", distribution))))[["q"]],
      method       = "L-BFGS-B",
      verbose      = verbose
    )

    if(abs(parameters[["par"]][1] - search_bounds1[1]) < .1 | abs(parameters[["par"]][1] - search_bounds1[2]) < .1)
      warning("The calibrated prior distribution is close to search boundary. Consider specifying a wider search range.")

    if(abs(parameters[["par"]][2] - search_bounds2[1]) < .1 | abs(parameters[["par"]][2] - search_bounds2[2]) < .1)
      warning("The calibrated prior distribution is close to search boundary. Consider specifying a wider search range.")

    parameters <- list(
      intercept = parameters[["par"]][1],
      aux       = parameters[["par"]][2]
    )

  }else{

    parameters <- optim(
      par   = 0,
      fn    = .solve_quantiles1,
      lower = search_bounds1[1],
      upper = search_bounds1[2],
      median_t     = median_t,
      quant_fun    = eval(parse(text = paste0(".", gsub("-", "_", distribution))))[["q"]],
      method       = "L-BFGS-B",
      verbose      = verbose
    )

    if(abs(parameters[["par"]][1] - search_bounds1[1]) < .1 | abs(parameters[["par"]][1] - search_bounds1[2]) < .1)
      warning("The calibrated prior distribution is close to search boundary. Consider specifying a wider search range.")

    parameters <- list(
      intercept = parameters[["par"]][1]
    )

  }

  if(verbose){
    cat("Parameters: \n")
    print(parameters)
    cat("\n\n")
  }

  return(parameters)
}

.solve_quantiles1 <- function(x, median_t, quant_fun, verbose){

  intercept <- x[1]

  res_median_t <- do.call(quant_fun, list(.5, eta = intercept))

  MSE <- (res_median_t - median_t)^2

  if(verbose){
    print(data.frame(intercept, res_median_t, MSE))
  }

  return(MSE)
}
.solve_quantiles2 <- function(x, median_t, iq_range_t, quant_fun, verbose){

  intercept <- x[1]
  aux       <- x[2]

  res_median_t   <- do.call(quant_fun, list(.50, eta = intercept, aux))
  res_25q_t      <- do.call(quant_fun, list(.25, eta = intercept, aux))
  res_75q_t      <- do.call(quant_fun, list(.75, eta = intercept, aux))
  res_iq_range_t <- res_75q_t - res_25q_t

  MSE <- (res_median_t - median_t)^2 + (res_iq_range_t - iq_range_t)^2

  if(verbose){
    print(data.frame(distribution, intercept, aux, res_median_t, res_iq_range_t, MSE))
  }

  if(is.na(MSE)){
    MSE <- 1e100
  }

  return(MSE)
}

.make_prior   <- function(parameters, prior_sd, distribution){

  if(.has_aux(distribution)){

    priors <- list(
      intercept = prior("normal",    list(parameters[["intercept"]], prior_sd[1])),
      aux       = prior("lognormal", list(lnorm_meanlog(parameters[["aux"]], prior_sd[2]), lnorm_sdlog(parameters[["aux"]], prior_sd[2])))
    )

  }else{

    priors <- list(
      intercept = prior("normal",    list(parameters[["intercept"]], prior_sd[1]))
    )

  }

  return(priors)
}

#' @title Create meta-analytic predictive prior distributions
#'
#' @description Calibrates prior distributions for parametric
#' survival analysis based on historical data.
#' Returns a list of prior distribution for the intercepts
#' and auxiliary parameters.
#'
#' @param datasets list of data.frames containing the historical
#' data. Each data.frame must contain a column named \code{"time"}
#' with the survival times and a column named \code{"status"}
#' with the censoring status.
#' @param distributions vector of parametric families for which
#' prior distributions ought to be calibrated
#' @param mu_scale scale of the prior Cauchy distribution for the meta-analytic
#' mean parameter
#' @param tau_scale scale of the prior Cauchy distribution for the meta-analytic
#' heterogeneity parameter
#'
#' @export
calibrate_meta_analytic  <- function(datasets, distributions = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"), mu_scale = 100, tau_scale = 10){

  if(!is.list(datasets))
    stop("'datasets' must be a list")
  if(!all(sapply(datasets, is.data.frame)))
    stop("'datasets' must be a list of data.frames")
  if(!all(sapply(datasets, function(df) all(c("time", "status") %in% colnames(df)))))
    stop("'datasets' must be a list of data.frames that contains 'time' and 'status' columns describing the survival times and censoring")

  # prepare output
  priors_int <- list()
  priors_aux <- list()

  # fit the survival models with flexsurv
  distributions <- unname(sapply(distributions, .RoBSA2flex_distribution))

  surv_fits <- list()
  for(df in seq_along(datasets)){
    surv_fits[[df]] <- list()
    for (d in distributions){
      surv_fits[[df]][[d]] <- suppressWarnings(flexsurv::flexsurvreg(Surv(time, status) ~ 1, data = datasets[[df]][datasets[[df]]$time != 0,], dist = d))
    }
  }

  coefs    <- sapply(distributions, function(fam)data.frame(do.call(rbind, sapply(surv_fits, function(df)coef(df[[fam]]), simplify = FALSE))), simplify = FALSE)
  coefs_se <- sapply(distributions, function(fam)data.frame(do.call(rbind, sapply(surv_fits, function(df)sqrt(diag(df[[fam]]$cov)), simplify = FALSE))), simplify = FALSE)
  colnames(coefs$exp)    <- "rate"
  colnames(coefs_se$exp) <- "rate"

  # fit meta-analytic priors
  coefs_meta <- list()
  for(fam in distributions){
    for(par in names(coefs[[fam]])){

      require("metaBMA")
      temp_fit <- suppressWarnings(metaBMA::meta_random(
        y   = coefs[[fam]][[par]], SE = coefs_se[[fam]][[par]],
        d   = metaBMA::prior("cauchy", c(location = 0, scale = mu_scale)),
        tau = metaBMA::prior("cauchy", c(location = 0, scale = tau_scale), lower = 0),
        control = list(adapt_delta = .95), cores = 1, chains = 2, iter = 5000, warmup = 2000))


      # deal with the negative parametrization for exp and gamma family
      if(fam %in% c("exp", "gamma") & par == "rate"){

        priors_int[[.flex2RoBSA_distribution(fam)]] <- prior(
          distribution = "normal",
          parameters   = list(
            mean = - temp_fit$estimates["d", "mean"],
            sd   = sqrt(temp_fit$estimates["d", "sd"]^2 + temp_fit$estimates["tau", "mean"]^2)
          ))

      }else{

        if(par %in% c("rate", "scale", "meanlog")){

          priors_int[[.flex2RoBSA_distribution(fam)]] <- prior(
            distribution = "normal",
            parameters   = list(
              mean = temp_fit$estimates["d", "mean"],
              sd   = sqrt(temp_fit$estimates["d", "sd"]^2 + temp_fit$estimates["tau", "mean"]^2)
          ))

        }else{

          priors_aux[[.flex2RoBSA_distribution(fam)]] <- prior(
            distribution = "lognormal",
            parameters   = list(
              meanlog = temp_fit$estimates["d", "mean"],
              sdlog   = sqrt(temp_fit$estimates["d", "sd"]^2 + temp_fit$estimates["tau", "mean"]^2)
          ))

        }

      }
    }
  }

  return(list(
    intercept = priors_int,
    aux       = priors_aux
  ))
}

.RoBSA2flex_distribution <- function(distribution){
  return(switch(
    distribution,
    "exp-aft"     = "exp",
    "weibull-aft" = "weibull",
    "lnorm-aft"   = "lognormal",
    "llogis-aft"  = "llogis",
    "gamma-aft"   = "gamma"
  ))
}
.flex2RoBSA_distribution <- function(distribution){
  return(switch(
    distribution,
    "exp"       = "exp-aft",
    "weibull"   = "weibull-aft",
    "lognormal" = "lnorm-aft",
    "llogis"    = "llogis-aft",
    "gamma"     = "gamma-aft"
  ))
}
