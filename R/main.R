#' @title Fit Robust Bayesian Survival Analysis
#'
#' @description Fits a RoBSA model. Please note
#' that the fitting function only supports dummy
#' coded predictors and cannot deal with factors.
#'
#' @param formula formula for the survival model
#' @param data data frame containing the data
#' @param distributions distributions of parametric
#' survival models
#' @param test_predictors vector of predictor names
#' to be tested with Bayesian model-averaged testing.
#' Defaults to \code{NULL}, no parameters are tested.
#' @param distributions_weights prior odds for the competing
#' distributions
#' @param prior_beta_null named list containing null prior
#' distribution for the predictors (with names corresponding
#' to the predictors)
#' @param prior_beta_alt named list containing prior
#' distribution for the predictors (with names corresponding
#' to the predictors)
#' @param prior_intercept named list containing prior
#' distribution for the intercepts (with names corresponding
#' to the distributions)
#' @param prior_aux named list containing prior
#' distribution for the auxiliary parameters (with names corresponding
#' to the distributions)
#' @param chains a number of chains of the MCMC algorithm.
#' @param sample a number of sampling iterations of the MCMC algorithm.
#' Defaults to \code{5000}.
#' @param burnin a number of burnin iterations of the MCMC algorithm.
#' Defaults to \code{2000}.
#' @param adapt a number of adaptation iterations of the MCMC algorithm.
#' Defaults to \code{500}.
#' @param thin a thinning of the chains of the MCMC algorithm. Defaults to
#' \code{1}.
#' @param parallel whether the individual models should be fitted in parallel.
#' Defaults to \code{FALSE}. The implementation is not completely stable
#' and might cause a connection error.
#' @param autofit whether the model should be fitted until the convergence
#' criteria (specified in \code{autofit_control}) are satisfied. Defaults to
#' \code{TRUE}.
#' @param autofit_control allows to pass autofit control settings with the
#' [set_autofit_control()] function. See \code{?set_autofit_control} for
#' options and default settings.
#' @param convergence_checks automatic convergence checks to assess the fitted
#' models, passed with [set_convergence_checks()] function. See
#' \code{?set_convergence_checks} for options and default settings.
#' @param save whether all models posterior distributions should be kept
#' after obtaining a model-averaged result. Defaults to \code{"all"} which
#' does not remove anything. Set to \code{"min"} to significantly reduce
#' the size of final object, however, some model diagnostics and further
#' manipulation with the object will not be possible.
#' @param seed a seed to be set before model fitting, marginal likelihood
#' computation, and posterior mixing for reproducibility of results. Defaults
#' to \code{NULL} - no seed is set.
#' @param silent whether all print messages regarding the fitting process
#' should be suppressed. Defaults to \code{TRUE}. Note that \code{parallel = TRUE}
#' also suppresses all messages.
#' @param ... additional arguments.
#'
#' @rdname RoBSA
#' @aliases RoBSA
#' @export
RoBSA <- function(
  formula, data, priors = NULL, test_predictors = NULL,

  distributions = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"),
  distributions_weights = rep(1, length(distributions)),

  # default prior distribution
  prior_beta_null   = get_default_prior_beta_null(),
  prior_beta_alt    = get_default_prior_beta_alt(),
  prior_factor_null = get_default_prior_factor_null(),
  prior_factor_alt  = get_default_prior_factor_alt(),
  prior_intercept   = get_default_prior_intercept(),
  prior_aux         = get_default_prior_aux(),

  # MCMC fitting settings
  chains = 3, sample = 5000, burnin = 2000, adapt = 500, thin = 1, parallel = FALSE,
  autofit = TRUE, autofit_control = set_autofit_control(), convergence_checks = set_convergence_checks(),

  # additional settings
  save = "all", seed = NULL, silent = TRUE, ...){

  dots         <- .RoBSA_collect_dots(...)
  object       <- NULL
  object$call  <- match.call()

  ### prepare & check the data
  object$data    <- .prepare_data(formula, data)
  object$formula <- formula

  ### check MCMC settings
  object$fit_control        <- BayesTools::JAGS_check_and_list_fit_settings(chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)
  object$autofit_control    <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = autofit_control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks)

  ### prepare and check the settings
  object$priors  <- .check_and_list_priors(priors = priors, distributions = distributions, data = object[["data"]], test_predictors = test_predictors,
                                           default_prior_beta_null = prior_beta_null, default_prior_beta_alt = prior_beta_alt,
                                           default_prior_factor_null = prior_factor_null, default_prior_factor_alt = prior_factor_alt,
                                           default_prior_intercept = prior_intercept, default_prior_aux = prior_aux)
  object$models  <- .prepare_models(object$priors, distributions, distributions_weights)


  ### fit the models and compute marginal likelihoods
  if(!object$fit_control[["parallel"]]){

    if(dots[["is_JASP"]]){
      .JASP_progress_bar_start(length(object[["models"]]))
    }

    for(i in seq_along(object[["models"]])){
      object$models[[i]] <- .fit_RoBSA_model(object, i)
      if(dots[["is_JASP"]]){
        .JASP_progress_bar_tick()
      }
    }

  }else{

    fitting_order <- .fitting_priority(object[["models"]])

    cl <- parallel::makePSOCKcluster(floor(RoBSA.get_option("max_cores") / object$fit_control[["chains"]]))
    parallel::clusterEvalQ(cl, {library("RoBSA")})
    parallel::clusterExport(cl, "object", envir = environment())
    object$models <- parallel::parLapplyLB(cl, fitting_order, .fit_RoBSA_model, object = object)[order(fitting_order)]
    parallel::stopCluster(cl)

  }

  return(object)
    # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)


  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){

    # TODO: any idea how to proceed with this?
    # balance probability of non-converged models TODO
    # if(object$control$balance_prob & any(!object$add_info$converged))object <- .balance_prob(object, object$add_info$converged)


    ### compute the model-space results
    object$RoBSA         <- .model_inference(object)
    object$coefficients  <- .compute_coeficients(object$RoBSA)
  }


  ### add warnings
  object$add_info$warnings <- c(object$add_info$warnings, .model_refit_warnings(sapply(1:length(object$models), function(i)object$models[[i]]$metadata, simplify = FALSE)))
  object$add_info$warnings <- c(object$add_info$warnings, .model_convergence_warnings(object))


  ### remove model posteriors if asked to
  if(save == "min"){
    for(i in 1:length(object$models)){
      if(length(object$models[[1]]$fit) != 1){
        object$models[[i]]$fit$mcmc <- NULL
      }
    }
  }


  ### print warnings
  if(!is.null(object$add_info$warnings)){
    for(w in object$add_info$warnings)
      warning(w)
  }
  if(sum(!object$add_info$converged) > 0)
    warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))


  class(object) <- "RoBSA"
  return(object)
}


update.RoBSA <- function(object, refit_failed = TRUE,
                         priors = NULL, prior_odds = NULL,
                         control = NULL, chains = NULL, iter = NULL, burnin = NULL, thin = NULL, parallel = NULL, seed = NULL, ...){


  # if(object$add_info$save == "min")
  #   stop("Models cannot be updated because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' in while fitting the model (see ?RoBSA for more details).")



  ### choose proper action based on the supplied input
  if(!is.null(priors)){

    what_to_do <- "fit_new_model"

    if(!is.null(prior_odds))object$models[[length(object$models)]]$prior_odds     <- prior_odds
    if(!is.null(prior_odds))object$models[[length(object$models)]]$prior_odds_set <- prior_odds


  }else if(!is.null(prior_odds)){

    what_to_do <- "update_prior_odds"
    if(length(prior_odds) != length(object$models))
      stop("The number of newly specified prior odds does not match the number of models. See '?update.RoSMA' for more details.")
    for(i in 1:length(object$models)){
      object$models[[i]]$prior_odds     <- prior_odds[i]
      object$models[[i]]$prior_odds_set <- prior_odds[i]
    }

  }else if(refit_failed & any(!object$add_info$converged)){

    what_to_do <- "refit_failed_models"

  }else{

    what_to_do <- "update_settings"

  }


  ### update control settings if any change is specified


  ### do the stuff
  if(what_to_do == "fit_new_model"){

    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    object$models[[length(object$models)]] <- .fit_RoBSA_wrap(object, length(object$models))

  }else if(what_to_do == "refit_failed_models"){

    converged_models <- .get_converged_models(object)
    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    for(i in c(1:length(object$models))[!converged_models]){
      object$models[[i]] <- .fit_RoBSA_wrap(object, i)
    }

  }else if(what_to_do == "transform_estimates"){

    for(i in c(1:length(object$models))){
      object$models[[i]] <- .transform_posteriors(object$models[[i]], object$add_info$output_scale, .transformation_var(output_scale))
    }
    object <- .transform_posterior(object, object$add_info$output_scale, .transformation_var(output_scale))
    object$add_info$output_scale <- .transformation_var(output_scale)

    return(object)

  }


  # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)

  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){

    ### compute the model-space results
    object$RoBSA         <- .model_inference(object)
    object$coefficients  <- .compute_coeficients(object$RoBSA)
  }


  ### add warnings
  object$add_info$warnings <- c(object$add_info$warnings, .model_refit_warnings(sapply(1:length(object$models), function(i)object$models[[i]]$metadata, simplify = FALSE)))
  object$add_info$warnings <- c(object$add_info$warnings, .model_convergence_warnings(object))


  ### remove model posteriors if asked to
  # if(save == "min"){
  #   for(i in 1:length(object$models)){
  #     if(length(object$models[[1]]$fit) != 1){
  #       object$models[[i]]$fit$mcmc <- NULL
  #     }
  #   }
  # }


  ### print warnings
  if(!is.null(object$add_info$warnings)){
    for(w in object$add_info$warnings)
      warning(w)
  }
  if(sum(!object$add_info$converged) > 0)
    warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))

  return(object)
}



.get_converged_models       <- function(object){

  converged <- NULL

  # basic convergence checks
  for(i in 1:length(object$models)){
    if(any(class(object$models[[i]]$fit) %in% c("simpleError", "error")) | is.infinite(object$models[[i]]$marg_lik$logml) | is.na(object$models[[i]]$marg_lik$logml)){
      converged <- c(converged, FALSE)
    }else{
      converged <- c(converged, TRUE)
    }
  }

  object$models <- object$models[converged]

  # remove models with unsatisfactory performance
  if(!is.null(object$control$allow_max_error) |!is.null(object$control$allow_max_rhat) | !is.null(object$control$allow_min_ESS)){
    diagnostics_summary <- summary.RoBSA(object, type = "models", diagnostics = TRUE, include_theta = object$control$allow_inc_theta)$diagnostics

    # deal with NAs for null models
    diagnostics_summary$"max(MCMC error)"[is.na(diagnostics_summary$"max(MCMC error)")] <- 0
    diagnostics_summary$"max(Rhat)"[is.na(diagnostics_summary$"max(Rhat)")]             <- 0
    diagnostics_summary$"min(ESS)"[is.na(diagnostics_summary$"min(ESS)")]               <- Inf


    if(!is.null(object$control$allow_max_error)){
      converged <- converged & (diagnostics_summary$"max(MCMC error)" < object$control$allow_max_error)
    }
    if(!is.null(object$control$allow_max_Rhat)){
      converged <- converged & diagnostics_summary$"max(Rhat)" < object$control$allow_max_rhat
    }
    if(!is.null(object$control$allow_min_ESS)){
      converged <- converged & diagnostics_summary$"min(ESS)"  > object$control$allow_min_ESS
    }
  }

  return(converged)
}
.model_refit_warnings       <- function(metadata){

  new_warn <- NULL

  # extract meta-data with fit-refit information
  refit_info <- t(sapply(metadata, function(x){
    if(is.null(x$refit_info)){
      return(c(x$i, NA))
    }else{
      return(c(x$i, x$refit_info))
    }
  }))

  marglik_info <- t(sapply(metadata, function(x){
    if(is.null(x$marg_lik)){
      return(c(x$i, NA))
    }else{
      return(c(x$i, x$marg_lik))
    }
  }))

  if(is.null(dim(refit_info)))  refit_info   <- matrix(refit_info,   ncol = 2)
  if(is.null(dim(marglik_info)))marglik_info <- matrix(marglik_info, ncol = 2)


  if(length(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Initial fit of %1$s %2$s failed due to incompatible starting values (most likely due to an outlier in the data and limited precision of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data.",
      ifelse(length(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1]) == 1, "model", "models"),
      paste(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Initial fit of %1$s %2$s failed due to incompatible starting values (most likely due to an outlier in the data and limited precision of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data and the model was refitted using boost likelihood function.",
      ifelse(length(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1]) == 1, "model", "models"),
      paste(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(refit_info[!refit_info[, 2] %in% c("empirical init", "refit with boost") & !is.na(refit_info[,2]), 1]) > 0){
    refit_info_messages_i <- refit_info[refit_info[, 2] != "not enough iterations" & !is.na(refit_info[,2]), 1]
    refit_info_messages   <- refit_info[refit_info[, 2] != "not enough iterations" & !is.na(refit_info[,2]), 2]

    for(i in 1:length(refit_info_messages_i)){
      new_warn <- c(new_warn, paste0("Model ", refit_info_messages_i[i]," failed with the following error: ", refit_info_messages[i]))
    }
  }


  if(length(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Marginal likelihood computation of %1$s %2$s couldn't be completed within the specified number of iterations.",
      ifelse(length(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1]) == 1, "model", "models"),
      paste(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 1]) > 0){
    marglik_info_messages_i <- marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 1]
    marglik_info_messages   <- marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 2]

    for(i in 1:length(marglik_info_messages_i)){
      new_warn <- c(new_warn, paste0("Marginal likelihood computation of model ", marglik_info_messages_i[i]," failed with the following error: ", marglik_info_messages[i]))
    }
  }

  return(new_warn)
}
.model_convergence_warnings <- function(object){

  new_warn <- NULL

  # used set values if specified by the user
  threshold_error <- ifelse(is.null(object$control$allow_max_error), Inf, object$control$allow_max_error)
  threshold_rhat  <- ifelse(is.null(object$control$allow_max_rhat), 1.05, object$control$allow_max_rhat)
  threshold_ESS   <- ifelse(is.null(object$control$allow_max_error), 100, object$control$allow_min_ESS)

  # get the diagnostics summary
  diagnostics_summary <- summary.RoBSA(object, type = "models", diagnostics = TRUE, include_theta = object$control$allow_inc_theta)$diagnostics

  # deal with NAs for null models
  diagnostics_summary$"max(MCMC error)"[is.na(diagnostics_summary$"max(MCMC error)")] <- 0
  diagnostics_summary$"max(Rhat)"[is.na(diagnostics_summary$"max(Rhat)")]             <- 0
  diagnostics_summary$"min(ESS)"[is.na(diagnostics_summary$"min(ESS)")]               <- Inf

  # find the problematic models
  warning_error <- rownames(diagnostics_summary)[diagnostics_summary$"max(MCMC error)" > threshold_error]
  warning_rhat  <- rownames(diagnostics_summary)[diagnostics_summary$"max(Rhat)"       > threshold_rhat]
  warning_ESS   <- rownames(diagnostics_summary)[diagnostics_summary$"min(ESS)"        < threshold_ESS]

  # add warnings messages
  if(length(warning_error) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with MCMC error larger than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_error) == 1, "Model", "Models"),
      paste(warning_error, collapse = ", "),
      threshold_error
    ))
  }

  if(length(warning_rhat) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with R-hat larger than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_rhat) == 1, "Model", "Models"),
      paste(warning_rhat, collapse = ", "),
      threshold_rhat
    ))
  }

  if(length(warning_ESS) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with ESS lower than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_ESS) == 1, "Model", "Models"),
      paste(warning_ESS, collapse = ", "),
      threshold_ESS
    ))
  }

  return(new_warn)
}
.inclusion_BF               <- function(prior_weights, posterior_weights, conditional_models){
  (sum(posterior_weights[conditional_models])/sum(posterior_weights[!conditional_models]))  /
    (sum(prior_weights[conditional_models])/sum(prior_weights[!conditional_models]))
}
.model_inference            <- function(object, n_samples = 10000){

  models    <- object$models
  data      <- object$data
  converged <- object$add_info$converged
  seed      <- object$control$seed

  predictors    <- attr(object$data, "predictors")
  distributions <- as.character(sapply(1:length(object$models), function(i)object$models[[i]]$distribution))

  # extract marginal likelihoods
  marg_liks <- sapply(models, function(x)x$marg_lik$logml)

  # determine the type of the models
  mm_predictors <- list()
  for(i in seq_along(predictors)){
    mm_predictors[[predictors[i]]] <- sapply(models, function(m)m$priors[["predictors"]][[predictors[i]]][["type"]] == "alt")
  }
  mm_distributions <- list()
  for(i in seq_along(unique(distributions))){
    mm_distributions[[unique(distributions)[i]]] <- sapply(models, function(m)m$distribution == unique(distributions)[i])
  }

  # extract model weights
  prior_weights_all        <- sapply(models, function(m)m$prior_odds)
  prior_weights_predictors <- list()
  for(i in seq_along(predictors)){
    prior_weights_predictors[[predictors[i]]] <- ifelse(mm_predictors[[predictors[i]]], prior_weights_all, 0)
  }
  prior_weights_distributions <- list()
  for(i in seq_along(unique(distributions))){
    prior_weights_distributions[[unique(distributions)[i]]] <- ifelse(mm_distributions[[unique(distributions)[i]]], prior_weights_all, 0)
  }

  # standardize model weights
  prior_weights_all   <- prior_weights_all   / sum(prior_weights_all)
  for(i in seq_along(predictors)){
    prior_weights_predictors[[predictors[i]]] <- prior_weights_predictors[[predictors[i]]]/sum(prior_weights_predictors[[predictors[i]]])
  }
  for(i in seq_along(unique(distributions))){
    prior_weights_distributions[[unique(distributions)[i]]] <- prior_weights_distributions[[unique(distributions)[i]]]/sum(prior_weights_distributions[[unique(distributions)[i]]])
  }

  ### compute model weights
  # overall
  weights_all        <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_all)
  weights_predictors <- list()
  for(i in seq_along(predictors)){
    if(any(mm_predictors[[predictors[i]]]) & all(!is.nan(prior_weights_predictors[[predictors[i]]]))){
      weights_predictors[[predictors[i]]] <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_predictors[[predictors[i]]])
    }
  }
  weights_distributions <- list()
  for(i in seq_along(unique(distributions))){
    if(any(mm_distributions[[unique(distributions)[i]]]) & all(!is.nan(prior_weights_distributions[[unique(distributions)[i]]]))){
      weights_distributions[[unique(distributions)[i]]] <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_distributions[[unique(distributions)[i]]])
    }
  }

  ### update the distributions and predictors trackers only to those present
  predictors    <- names(weights_predictors)
  distributions <- names(weights_distributions)

  ### compute inclusion BFs
  BF_predictors    <- list()
  for(i in seq_along(predictors)){
    BF_predictors[[predictors[i]]] <- .inclusion_BF(prior_weights_all, weights_all, mm_predictors[[predictors[i]]])
  }
  BF_distributions <- list()
  for(i in seq_along(unique(distributions))){
    BF_distributions[[unique(distributions)[i]]] <- .inclusion_BF(prior_weights_all, weights_all, mm_distributions[[unique(distributions)[i]]])
  }


  ### sample and mix the individual posteriors
  if(!is.null(seed))set.seed(seed)
  samples <- list()
  samples$predictors <- list()
  for(i in seq_along(predictors)){
    samples$predictors$averaged[[predictors[i]]]    <- .mix_samples(models, weights_all, converged, predictors[i], n_samples, seed)
  }
  for(i in seq_along(predictors)){
    samples$predictors$conditional[[predictors[i]]] <- .mix_samples(models, weights_predictors[[predictors[i]]], converged, predictors[i], n_samples, seed)
  }
  for(i in seq_along(unique(distributions))){
    samples$intercept[[unique(distributions)[i]]] <- .mix_samples(models, weights_distributions[[unique(distributions)[i]]], converged, "intercept", n_samples, seed)
  }
  for(i in seq_along(unique(distributions))){
    if(.has_aux(unique(distributions)[i])){
      samples$aux[[unique(distributions)[i]]] <- .mix_samples(models, weights_distributions[[unique(distributions)[i]]], converged, "aux", n_samples, seed)
    }
  }


  ### edit names
  names(weights_all)   <- names(models)


  ### compute prior and posterior probabilities for model types
  prior_prob_all            <- prior_weights_all
  prior_prob_distributions  <- sapply(unique(distributions), function(distribution)sum(prior_weights_all[mm_distributions[[distribution]]]))
  prior_prob_predictors     <- sapply(predictors,            function(predictor)   sum(prior_weights_all[mm_predictors[[predictor]]]))

  posterior_prob_all            <- weights_all
  posterior_prob_distributions  <- sapply(unique(distributions), function(distribution)sum(weights_all[mm_distributions[[distribution]]]))
  posterior_prob_predictors     <- sapply(predictors,            function(predictor)   sum(weights_all[mm_predictors[[predictor]]]))


  # return the results
  output <- list(
    samples        = samples,
    BF             = list(
      distributions = BF_distributions,
      predictors    = BF_predictors
    ),
    prior_prob     = list(
      all               = prior_prob_all,
      distributions     = prior_prob_distributions,
      predictors        = prior_prob_predictors
    ),
    posterior_prob = list(
      all               = posterior_prob_all,
      distributions     = posterior_prob_distributions,
      predictors        = posterior_prob_predictors
    )
  )
  return(output)
}
.mix_samples                <- function(models, weights, converged, parameter, n_samples, seed){

  is.predictor <- !parameter %in% c("intercept", "aux")

  if(!is.null(seed)) set.seed(seed) else set.seed(1)
  samples <- NULL

  for(i in c(1:length(models))[converged]){

    model_samples <- suppressWarnings(coda::as.mcmc(models[[i]]$fit))
    # deal with the possibility of intercept only exponential model
    if(!is.matrix(model_samples)){
      model_samples <- matrix(model_samples, ncol = 1)
      colnames(model_samples) <- "intercept"
    }


    # creates indexing - the set seet at the beggining makes sure that the samples for different predictors are correlated
    ind <- sample(nrow(model_samples), round(n_samples * weights[i]), replace = TRUE)
    if(length(ind) == 0)next

    if(is.predictor){
      samples <- c(samples,
                   if(models[[i]]$priors[["predictors"]][[parameter]]$distribution == "point"){
                     rep(models[[i]]$priors[["predictors"]][[parameter]]$parameters$location, round(n_samples * weights[i]))
                   }else{
                     model_samples[ind, paste0("beta_", parameter)]
                   })
    }else{
      samples <- c(samples,
                   if(models[[i]]$priors[[parameter]]$distribution == "point"){
                     rep(models[[i]]$priors[[parameter]]$parameters$location, round(n_samples * weights[i]))
                   }else{
                     model_samples[ind, parameter]
                   })
    }

  }


  return(samples)
}
.compute_coeficients        <- function(RoBSA){
  sapply(RoBSA$samples$predictors$averaged, mean)
}
.BF_format                  <- function(BF, BF01 = FALSE, logBF = FALSE){
  BF[is.nan(BF)] <- NA
  if(BF01){
    BF <- 1/BF
  }else{
    BF <- BF
  }
  if(logBF){
    BF <- log(BF)
  }
  return(BF)
}
