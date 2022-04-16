#' @title Prints summary of \code{"RoBSA"} ensemble implied by the specified priors
#'
#' @description \code{check_setup} prints summary of \code{"RoBSA"} ensemble
#' implied by the specified prior distributions. It is useful for checking
#' the ensemble configuration prior to fitting all of the models.
#'
#' @inheritParams RoBSA
#' @param models should the models' details be printed.
#' @param silent do not print the results.
#'
#' @return \code{check_setup} invisibly returns list of summary tables.
#'
#' @seealso [RoBSA()]
#' @export
check_setup <- function(
    formula, data, priors = NULL, test_predictors = NULL,
    distributions = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"),
    distributions_weights = rep(1, length(distributions)),
    prior_beta_null  = get_default_prior_beta_null(),
    prior_beta_alt   = get_default_prior_beta_alt(),
    prior_intercept  = get_default_prior_intercept(),
    prior_aux        = get_default_prior_aux(),

    # MCMC fitting settings
    chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
    autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

    # additional settings
    save = "all", seed = NULL, silent = TRUE, ...){

  # .....

  object <- list()
  object$priors   <- .check_and_list_priors(tolower(model_type), priors_effect_null, priors_effect, priors_heterogeneity_null, priors_heterogeneity, priors_bias_null, priors_bias, priors_rho, priors_rho_null, object$add_info[["prior_scale"]])
  object$models   <- .make_models(object[["priors"]], multivariate = FALSE)


  ### model types overview
  effect         <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  bias           <- sapply(object$models, function(model)!.is_component_null(model[["priors"]], "bias"))

  # obtain the parameter types
  weightfunctions <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.weightfunction)))
  PET             <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
  PEESE           <- sapply(object$models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

  # number of model types
  n_models    <- c(
    mu    = sum(effect),
    tau   = sum(heterogeneity),
    omega = sum(bias)
  )

  # extract model weights
  prior_weights   <- sapply(object$models, function(m)m$prior_weights)
  # standardize model weights
  prior_weights   <- prior_weights / sum(prior_weights)
  # conditional model weights
  models_prior <- c(
    mu    <- sum(prior_weights[effect]),
    tau   <- sum(prior_weights[heterogeneity]),
    omega <- sum(prior_weights[bias])
  )

  # create overview table
  components <- data.frame(
    "models"     = n_models,
    "prior_prob" = models_prior
  )
  rownames(components) <- c("Effect", "Heterogeneity", "Bias")

  class(components)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(components))
  attr(components, "type")      <- c("n_models", "prior_prob")
  attr(components, "rownames")  <- TRUE
  attr(components, "n_models")  <- length(object$models)
  attr(components, "title")     <- "Components summary:"
  attr(components, "footnotes") <- NULL
  attr(components, "warnings")  <- NULL

  object$components <- components

  ### model details
  if(models){
    priors_effect        <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$mu, silent = TRUE))
    priors_heterogeneity <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors$tau, silent = TRUE))
    priors_bias          <- sapply(1:length(object$models), function(i){
      if(weightfunctions[i]){
        print(object$models[[i]]$priors$omega, silent = TRUE)
      }else if(PET[i]){
        print(object$models[[i]]$priors$PET, silent = TRUE)
      }else if(PEESE[i]){
        print(object$models[[i]]$priors$PEESE, silent = TRUE)
      }else{
        ""
      }
    })
    prior_weights  <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_weights)
    prior_prob     <- prior_weights / sum(prior_weights)

    summary <- cbind.data.frame(
      "Model"         = 1:length(object$models),
      "Effect"        = priors_effect,
      "Heterogeneity" = priors_heterogeneity,
      "Bias"          = priors_bias,
      "prior_prob"    = prior_prob
    )
    class(summary)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(summary))
    attr(summary, "type")      <- c("integer", rep("prior", 3), "prior_prob")
    attr(summary, "rownames")  <- FALSE
    attr(summary, "title")     <- "Models overview:"
    attr(summary, "footnotes") <- NULL
    attr(summary, "warnings")  <- NULL

    object$summary <- summary
  }


  if(!silent){
    cat("Robust Bayesian meta-analysis (set-up)\n")
    print(components, quote = FALSE, right = TRUE)

    if(models){
      cat("\n")
      print(summary, quote = FALSE, right = TRUE)
    }
  }

  return(invisible(object))
}



#' @title Control MCMC fitting process
#'
#' @description Controls settings for the autofit
#' process of the MCMC JAGS sampler (specifies termination
#' criteria), and values for the convergence checks.
#'
#' @param max_Rhat maximum value of the R-hat diagnostic.
#' Defaults to \code{1.05}.
#' @param min_ESS minimum estimated sample size.
#' Defaults to \code{500}.
#' @param max_error maximum value of the MCMC error.
#' Defaults to \code{NULL}. Be aware that PEESE publication bias
#' adjustment can have estimates on different scale than the rest of
#' the output, resulting in relatively large max MCMC error.
#' @param max_SD_error maximum value of the proportion of MCMC error
#' of the estimated SD of the parameter.
#' Defaults to \code{NULL}.
#' @param max_time list with the time and unit specifying the maximum
#' autofitting process per model. Passed to \link[base]{difftime} function
#' (possible units are \code{"secs", "mins", "hours", "days", "weeks", "years"}).
#' Defaults to \code{list(time = 60, unit = "mins")}.
#' @param sample_extend number of samples to extend the fitting process if
#' the criteria are not satisfied.
#' Defaults to \code{1000}.
#' @param remove_failed whether models not satisfying the convergence checks should
#' be removed from the inference. Defaults to \code{FALSE} - only a warning is raised.
#' @param balance_probability whether prior model probability should be balanced
#' across the combinations of models with the same H0/H1 for effect / heterogeneity / bias
#' in the case of non-convergence. Defaults to \code{TRUE}.
#'
#'
#' @return \code{set_autofit_control} returns a list of autofit control settings
#' and \code{set_convergence_checks} returns a list of convergence checks settings.
#'
#' @export set_autofit_control
#' @export set_convergence_checks
#' @name RoBSA_control
#' @aliases set_autofit_control, set_convergence_checks
#'
#' @seealso [RoBSA], [update.RoBSA]
NULL

#' @rdname RoBSA_control
set_autofit_control     <- function(max_Rhat = 1.05, min_ESS = 500, max_error = NULL, max_SD_error = NULL, max_time = list(time = 60, unit = "mins"), sample_extend = 1000){

  autofit_settings <- list(
    max_Rhat      = max_Rhat,
    min_ESS       = min_ESS,
    max_error     = max_error,
    max_SD_error  = max_SD_error,
    max_time      = max_time,
    sample_extend = sample_extend
  )
  autofit_settings <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_settings, call = "Checking 'autofit_control':\n\t")

  return(autofit_settings)
}
#' @rdname RoBSA_control
set_convergence_checks  <- function(max_Rhat = 1.05, min_ESS = 500, max_error = NULL, max_SD_error = NULL, remove_failed = FALSE, balance_probability = TRUE){

  convergence_checks <- list(
    max_Rhat            = max_Rhat,
    min_ESS             = min_ESS,
    max_error           = max_error,
    max_SD_error        = max_SD_error,
    remove_failed       = remove_failed,
    balance_probability = balance_probability
  )
  # allows NULL arguments so it can be used in this way too
  convergence_checks <- .check_and_list_convergence_checks(convergence_checks)

  return(convergence_checks)
}


.update_fit_control     <- function(old_fit_control, chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed){

  if(is.null(chains)){
    chains <- old_fit_control[["chains"]]
  }
  if(is.null(adapt)){
    adapt  <- old_fit_control[["adapt"]]
  }
  if(is.null(burnin)){
    burnin <- old_fit_control[["burnin"]]
  }
  if(is.null(sample)){
    sample <- old_fit_control[["sample"]]
  }
  if(is.null(thin)){
    thin  <- old_fit_control[["thin"]]
  }
  if(is.null(autofit)){
    autofit  <- old_fit_control[["autofit"]]
  }
  if(is.null(parallel)){
    parallel <- old_fit_control[["parallel"]]
  }
  if(is.null(silent)){
    silent <- old_fit_control[["silent"]]
  }
  if(is.null(seed)){
    seed   <- old_fit_control[["seed"]]
  }

  new_fit_control <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)

  return(new_fit_control)
}
.update_autofit_control <- function(old_autofit_control, autofit_control){

  if(!is.null(autofit_control[["max_Rhat"]])){
    max_Rhat <- autofit_control[["max_Rhat"]]
  }else{
    max_Rhat <- old_autofit_control[["max_Rhat"]]
  }
  if(!is.null(autofit_control[["min_ESS"]])){
    min_ESS <- autofit_control[["min_ESS"]]
  }else{
    min_ESS <- old_autofit_control[["min_ESS"]]
  }
  if(!is.null(autofit_control[["max_error"]])){
    max_error <- autofit_control[["max_error"]]
  }else{
    max_error <- old_autofit_control[["max_error"]]
  }
  if(!is.null(autofit_control[["max_SD_error"]])){
    max_SD_error <- autofit_control[["max_SD_error"]]
  }else{
    max_SD_error <- old_autofit_control[["max_SD_error"]]
  }
  if(!is.null(autofit_control[["max_time"]])){
    max_time <- autofit_control[["max_time"]]
  }else{
    max_time <- old_autofit_control[["max_time"]]
  }
  if(!is.null(autofit_control[["sample_extend"]])){
    sample_extend <- autofit_control[["sample_extend"]]
  }else{
    sample_extend <- old_autofit_control[["sample_extend"]]
  }

  new_autofit_control <- set_autofit_control(max_Rhat = max_Rhat, min_ESS = min_ESS, max_error = max_error, max_SD_error = max_SD_error, max_time = max_time, sample_extend = sample_extend)
  new_autofit_control <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = new_autofit_control)

  return(new_autofit_control)
}
.update_convergence_checks <- function(old_convergence_checks, convergence_checks){

  if(!is.null(convergence_checks[["max_Rhat"]])){
    max_Rhat <- convergence_checks[["max_Rhat"]]
  }else{
    max_Rhat <- old_convergence_checks[["max_Rhat"]]
  }
  if(!is.null(convergence_checks[["min_ESS"]])){
    min_ESS <- convergence_checks[["min_ESS"]]
  }else{
    min_ESS <- old_convergence_checks[["min_ESS"]]
  }
  if(!is.null(convergence_checks[["max_error"]])){
    max_error <- convergence_checks[["max_error"]]
  }else{
    max_error <- old_convergence_checks[["max_error"]]
  }
  if(!is.null(convergence_checks[["max_SD_error"]])){
    max_SD_error <- convergence_checks[["max_SD_error"]]
  }else{
    max_SD_error <- old_convergence_checks[["max_SD_error"]]
  }
  if(!is.null(convergence_checks[["remove_failed"]])){
    remove_failed <- convergence_checks[["remove_failed"]]
  }else{
    remove_failed <- old_convergence_checks[["remove_failed"]]
  }
  if(!is.null(convergence_checks[["balance_probability"]])){
    balance_probability <- convergence_checks[["balance_probability"]]
  }else{
    balance_probability <- old_convergence_checks[["balance_probability"]]
  }

  new_convergence_checks <- set_convergence_checks(max_Rhat = max_Rhat, min_ESS = min_ESS, max_error = max_error, max_SD_error = max_SD_error, remove_failed = remove_failed, balance_probability = balance_probability)
  new_convergence_checks <- .check_and_list_convergence_checks(new_convergence_checks)
}


.check_and_list_convergence_checks <- function(convergence_checks){

  remove_failed       <- convergence_checks[["remove_failed"]]
  balance_probability <- convergence_checks[["balance_probability"]]
  convergence_checks["remove_failed"]       <- NULL
  convergence_checks["balance_probability"] <- NULL
  convergence_checks <- BayesTools::JAGS_check_and_list_autofit_settings(convergence_checks, skip_sample_extend = TRUE, call = "Checking 'convergence_checks':\n\t")

  BayesTools::check_bool(remove_failed,       "remove_failed",       call = "Checking 'convergence_checks':\n\t")
  BayesTools::check_bool(balance_probability, "balance_probability", call = "Checking 'convergence_checks':\n\t")
  convergence_checks[["remove_failed"]]       <- remove_failed
  convergence_checks[["balance_probability"]] <- balance_probability
  return(convergence_checks)
}
.check_and_list_add_info           <- function(distributions, predictors, predictors_test, seed, save, rescale_data, warnings, errors){

  BayesTools::check_char(distributions, "distributions", allow_values = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"), check_length = FALSE)
  BayesTools::check_char(predictors, "predictors", allow_NULL = TRUE, check_length = FALSE)
  BayesTools::check_char(predictors_test, "predictors_test", allow_NULL = TRUE, check_length = FALSE)
  BayesTools::check_real(seed, "seed", allow_NULL = TRUE)
  BayesTools::check_char(save, "save", allow_values = c("min", "all"))
  BayesTools::check_bool(rescale_data, "rescale_data")

  return(list(
    distributions    = distributions,
    predictors       = predictors,
    predictors_test  = predictors_test,
    seed             = seed,
    save             = save,
    rescale_data     = rescale_data,
    warnings         = warnings,
    errors           = errors
  ))
}

#' @title Default prior distributions
#'
#' @description Functions for setting default prior distributions.
#'
#' @return \code{get_default_prior_beta_null} and \code{get_default_prior_beta_alt}
#' return a prior distribution and \code{get_default_prior_intercept} and
#' \code{get_default_prior_aux} return a list of prior distributions.
#'
#' @export get_default_prior_beta_null
#' @export get_default_prior_beta_alt
#' @export get_default_prior_factor_null
#' @export get_default_prior_factor_null
#' @export get_default_prior_intercept
#' @export get_default_prior_aux
#' @name default_prior
#' @aliases get_default_prior_beta_null, get_default_prior_beta_alt,
#' get_default_prior_factor_null, get_default_prior_factor_alt,
#' get_default_prior_intercept, get_default_prior_aux
NULL

#' @rdname default_prior
get_default_prior_beta_null <- function(){
  prior("spike", list(location = 0))
}
#' @rdname default_prior
get_default_prior_beta_alt  <- function(){
  prior("normal", list(mean = 0, sd = 1))
}
#' @rdname default_prior
get_default_prior_factor_null  <- function(){
  list(
    "treatment"   = prior("spike", parameters = list("location" = 0)),
    "orthonormal" = prior("spike", parameters = list("location" = 0))
  )
}
#' @rdname default_prior
get_default_prior_factor_alt   <- function(){
  list(
    "treatment"   = prior_factor("normal",  list(mean = 0, sd = 1), contrast = "treatment"),
    "orthonormal" = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
}
#' @rdname default_prior
get_default_prior_intercept <- function(){
  list(
    "exp-aft"     = prior("normal", list(mean = 0, sd = 5)),
    "weibull-aft" = prior("normal", list(mean = 0, sd = 5)),
    "lnorm-aft"   = prior("normal", list(mean = 0, sd = 5)),
    "llogis-aft"  = prior("normal", list(mean = 0, sd = 5)),
    "gamma-aft"   = prior("normal", list(mean = 0, sd = 5))
  )
}
#' @rdname default_prior
get_default_prior_aux       <- function(){
  list(
    "exp-aft"     = NULL,
    "weibull-aft" = prior("normal", list(mean = 0, sd = 1), list(0, Inf)),
    "lnorm-aft"   = prior("normal", list(mean = 0, sd = 1), list(0, Inf)),
    "llogis-aft"  = prior("normal", list(mean = 0, sd = 1), list(0, Inf)),
    "gamma-aft"   = prior("normal", list(mean = 0, sd = 1), list(0, Inf))
  )
}

# some functions for the JASP implementation
.RoBSA_collect_dots      <- function(...){

  dots <- list(...)

  known_dots <- c("is_JASP")
  if(any(!names(dots) %in% known_dots))
    stop(paste0("Uknown arguments to 'RoBSA': ", paste("'", names(dots)[!names(dots) %in% known_dots], "'", collapse = ", "), "."), call. = FALSE)

  if(is.null(dots[["is_JASP"]])){
    dots[["is_JASP"]] <- FALSE
  }else{
    dots[["is_JASP"]] <- TRUE
  }

  return(dots)
}
.JASP_progress_bar_start <- function(n){
  eval(expr = parse(text = 'startProgressbar(n)'))
}
.JASP_progress_bar_tick  <- function(){
  eval(expr = parse(text = 'progressbarTick()'))
}
