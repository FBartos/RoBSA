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
#' @param rescale_data whether continuous predictors should be rescaled prior to
#' estimating the model. Defaults to \code{FALSE}.
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
  save = "all", seed = NULL, silent = TRUE, rescale_data = FALSE, ...){

  dots         <- .RoBSA_collect_dots(...)
  object       <- NULL
  object$call  <- match.call()


  ### prepare & check the data
  object$data    <- .prepare_data(formula, data, rescale_data)
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


  ### additional information
  object$add_info <- .check_and_list_add_info(
    distributions    = distributions,
    predictors       = attr(object[["priors"]], "terms"),
    predictors_test  = attr(object[["priors"]], "terms_test"),
    seed             = seed,
    save             = save,
    rescale_data     = rescale_data,
    warnings         = attr(object[["data"]], "warnings"),
    errors           = NULL
  )


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


  # create ensemble only if at least one model converged
  if(any(.get_model_convergence(object))){

    # TODO? balance probability of non-converged models?

    ### compute the model-space results
    object$models        <- BayesTools::models_inference(object[["models"]])
    object$RoBSA         <- .ensemble_inference(object)
    object$coefficients  <- .compute_coeficients(object[["RoBSA"]])
  }


  ### collect and print errors and warnings
  object$add_info[["errors"]]   <- c(object$add_info[["errors"]],   .get_model_errors(object))
  object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], .get_model_warnings(object))
  .print_errors_and_warnings(object)


  ### remove model posteriors if asked to
  if(save == "min"){
    object <- .remove_model_posteriors(object)
    object <- .remove_model_margliks(object)
  }


  class(object) <- "RoBSA"
  return(object)
}



#' @title Updates a fitted RoBSA object
#'
#' @description \code{update.RoBSA} can be used to
#' \enumerate{
#'   \item{add an additional model to an existing \code{"RoBSA"} object by
#'    specifying the distribution, and either null or alternative priors
#'    for each parameter and prior weight of the model,}
#'   \item{change the prior weights of fitted models by specifying a vector
#'   \code{prior_weights} of the same length as the fitted models,}
#'   \item{refitting models that failed to converge with updated settings
#'   of control parameters,}
#'   \item{or changing the convergence criteria and recalculating the ensemble
#'   results by specifying new \code{control} argument and setting
#'   \code{refit_failed == FALSE}.}
#' }
#'
#' @param object a fitted RoBSA object
#' @param distribution a distribution of the new model.
#' @param model_weights either a single value specifying prior model weight
#' of a newly specified model using priors argument, or a vector of the
#' same length as already fitted models to update their prior weights.
#' @param refit_failed whether failed models should be refitted. Relevant only
#' if new priors or \code{prior_weights} are not supplied. Defaults to \code{TRUE}.
#' @inheritParams RoBSA
#' @param ... additional arguments.
#'
#' @details See [RoBSA()] for more details.
#'
#'
#'
#' @return \code{RoBSA} returns an object of class 'RoBSA'.
#'
#' @seealso [RoBSA()], [summary.RoBSA()], [prior()], [check_setup()]
#' @export
update.RoBSA <- function(object, refit_failed = TRUE,

                         distribution  = NULL,
                         model_weights = NULL,

                         # default prior distribution
                         prior_beta_null   = get_default_prior_beta_null(),
                         prior_beta_alt    = get_default_prior_beta_alt(),
                         prior_factor_null = get_default_prior_factor_null(),
                         prior_factor_alt  = get_default_prior_factor_alt(),
                         prior_intercept   = get_default_prior_intercept(),
                         prior_aux         = get_default_prior_aux(),

                         chains = NULL, adapt = NULL, burnin = NULL, sample = NULL, thin = NULL, autofit = NULL, parallel = NULL,
                         autofit_control = NULL, convergence_checks = NULL,
                         save = "all", seed = NULL, silent = TRUE, ...){

  # TODO: add ability to change the output scale
  output_scale <- NULL

  if(object$add_info$save == "min")
    stop("Models cannot be updated because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' in while fitting the model (see ?RoBSA for more details).")


  ### choose proper action based on the supplied input
  if(!is.null(distribution)){

    what_to_do <- "fit_new_model"
    new_priors  <- .check_and_list_priors(priors = priors, distributions = distributions, data = object[["data"]], test_predictors = test_predictors,
                                          default_prior_beta_null = prior_beta_null, default_prior_beta_alt = prior_beta_alt,
                                          default_prior_factor_null = prior_factor_null, default_prior_factor_alt = prior_factor_alt,
                                          default_prior_intercept = prior_intercept, default_prior_aux = prior_aux)
    object$models[length(object$models) + 1]  <- list(.prepare_models(new_priors, distributions, distributions_weights)[[1]])

    if(!is.null(model_weights)){
      object$models[[length(object$models)]]$model_weights     <- model_weights
      object$models[[length(object$models)]]$model_weights_set <- model_weights
    }


  }else if(!is.null(model_weights)){

    what_to_do <- "update_model_weights"
    if(length(model_weights) != length(object$models))
      stop("The number of newly specified prior odds does not match the number of models. See '?update.RoBSA' for more details.")
    for(i in 1:length(object$models)){
      object$models[[i]]$prior_weights     <- model_weights[i]
      object$models[[i]]$prior_weights_set <- model_weights[i]
    }

  }else if(refit_failed & any(!.get_model_convergence(object))){

    what_to_do <- "refit_failed_models"

  }else{

    what_to_do <- "update_settings"

  }


  ### update control settings if any change is specified
  object[["fit_control"]]        <- .update_fit_control(object[["fit_control"]], chains = chains, adapt = adapt, burnin = burnin, sample = sample, thin = thin, autofit = autofit, parallel = parallel, cores = chains, silent = silent, seed = seed)
  object[["autofit_control"]]    <- .update_autofit_control(object[["autofit_control"]], autofit_control)
  object[["convergence_checks"]] <- .update_convergence_checks(object[["convergence_checks"]], convergence_checks)


  ### do the stuff
  if(what_to_do == "fit_new_model"){

    object[["models"]][[length(object$models)]] <- .fit_RoBSA_model(object, length(object$models))

  }else if(what_to_do == "refit_failed_models"){

    for(i in c(1:length(object$models))[!.get_model_convergence(object)]){
      object[["models"]][[i]] <- .fit_RoBSA_model(object, i)
    }

  }


  # create ensemble only if at least one model converged
  if(any(.get_model_convergence(object))){

    # TODO? balance probability of non-converged models?

    ### compute the model-space results
    object$models        <- BayesTools::models_inference(object[["models"]])
    object$RoBSA         <- .ensemble_inference(object)
    object$coefficients  <- .compute_coeficients(object[["RoBSA"]])
  }


  ### collect and print errors and warnings
  object$add_info[["errors"]]   <- c(object$add_info[["errors"]],   .get_model_errors(object))
  object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], .get_model_warnings(object))
  .print_errors_and_warnings(object)


  ### remove model posteriors if asked to
  if(save == "min"){
    object <- .remove_model_posteriors(object)
    object <- .remove_model_margliks(object)
  }

  return(object)
}
