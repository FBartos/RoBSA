# the main functions
.fit_RoBSA_model <- function(object, i){

  model              <- object[["models"]][[i]]
  priors             <- model[["priors"]]
  fit_control        <- object[["fit_control"]]
  autofit_control    <- object[["autofit_control"]]
  convergence_checks <- object[["convergence_checks"]]
  formula            <- object[["formula"]]
  data               <- object[["data"]]

  errors   <- NULL
  warnings <- NULL

  if(!fit_control[["silent"]]){
    cat(paste0("\nFitting model [", i, "]\n"))
  }

  # don't sample the complete null model
  if(!.is_model_constant(priors)){

    # generate the model syntax
    model_syntax <- .generate_model_syntax(model[["distribution"]], data)

    formula_list        <- .generate_model_formula_list(formula)
    formula_data_list   <- .generate_model_formula_data_list(data)
    formula_prior_list  <- .generate_model_formula_prior_list(priors)

    # remove unnecessary objects from data to mitigate warnings
    fit_data   <- .fit_data(data)
    fit_priors <- .fit_priors(priors)

    # fit the model
    fit <- BayesTools::JAGS_fit(
      model_syntax       = model_syntax,
      data               = fit_data,
      prior_list         = fit_priors,
      formula_list       = formula_list,
      formula_data_list  = formula_data_list,
      formula_prior_list = formula_prior_list,
      chains             = fit_control[["chains"]],
      adapt              = fit_control[["adapt"]],
      burnin             = fit_control[["burnin"]],
      sample             = fit_control[["sample"]],
      thin               = fit_control[["thin"]],
      autofit            = fit_control[["autofit"]],
      autofit_control    = autofit_control,
      parallel           = fit_control[["parallel"]],
      cores              = fit_control[["cores"]],
      silent             = fit_control[["silent"]],
      seed               = fit_control[["seed"]],
      required_packages  = "RoBSA"
    )

    # assess the model fit and deal with errors
    if(inherits(fit, "error")){

      if(grepl("Unknown function", fit$message))
        stop("The RoBSA module could not be loaded. Check whether the RoBMA package was installed correctly and whether 'RoBMA::RoBMA.private$module_location' contains path to the RoBMA JAGS module.")

      fit            <- list()
      converged      <- FALSE
      has_posterior  <- FALSE
      errors         <- c(errors, fit$message)
      # deal with failed models
      marglik        <- list()
      marglik$logml  <- NA
      class(marglik) <- "bridge"

    }else{

      has_posterior <- TRUE
      check_fit     <- BayesTools::JAGS_check_convergence(
        fit          = fit,
        prior_list   = attr(fit, "prior_list"),
        max_Rhat     = convergence_checks[["max_Rhat"]],
        min_ESS      = convergence_checks[["min_ESS"]],
        max_error    = convergence_checks[["max_error"]],
        max_SD_error = convergence_checks[["max_SD_error"]]
      )
      warnings    <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
      if(convergence_checks[["remove_failed"]] && !check_fit){
        converged <- FALSE
      }else{
        converged <- TRUE
      }

    }

    # compute marginal likelihood
    if(length(fit) != 0){

      marglik <- BayesTools::JAGS_bridgesampling(
        fit                = fit,
        log_posterior      = .marglik_function,
        data               = fit_data,
        prior_list         = fit_priors,
        formula_list       = formula_list,
        formula_data_list  = formula_data_list,
        formula_prior_list = formula_prior_list,
        maxiter            = 50000,
        silent             = fit_control[["silent"]],
        distribution       = model[["distribution"]]
      )

      # deal with failed marginal likelihoods
      if(inherits(marglik, "error")){

        errors         <- c(errors, marglik$message)
        converged      <- FALSE
        marglik        <- list()
        marglik$logml  <- NA
        class(marglik) <- "bridge"

      }else{

        # forward warnings if present
        warnings <- c(warnings, attr(marglik, "warnings"))

      }
    }


  }else{

    fit_data                <- .fit_data(object[["data"]], priors, add_info[["effect_direction"]], add_info[["prior_scale"]])
    converged               <- TRUE
    has_posterior           <- FALSE
    fit                     <- list()
    attr(fit, "prior_list") <- priors
    class(fit)              <- "null_model"
    marglik                 <- list()
    marglik$logml           <- .marglik_function_null(priors, data)
    class(marglik)          <- "bridge"

  }

  # add model summaries
  if(has_posterior){
    fit_summary <- BayesTools::runjags_estimates_table(fit = fit, warnings = warnings)
  }else{
    fit_summary <- BayesTools::runjags_estimates_empty_table()
  }

  model <- c(
    model,
    fit           = list(fit),
    fit_summary   = list(fit_summary),
    marglik       = list(marglik),
    errors        = list(errors),
    warnings      = list(warnings),
    converged     = converged,
    has_posterior = has_posterior
  )

  return(model)
}

# model parameter tools
.fit_data    <- function(data){

  types <- attr(data[["survival"]], "type")

  fit_data <- list("time" = unname(do.call(c, data[["survival"]][c(
    if(any(types == "event"))  "t_event",
    if(any(types == "cens_r")) "t_cens_r",
    if(any(types == "cens_l")) "t_cens_l",
    if(any(types == "cens_i")) c("t_cens_ir", "t_cens_il")
  )])))

  return(fit_data)
}
.fit_priors  <- function(priors){
  if(is.null(priors[["aux"]])){
    return(NULL)
  }else{
    return(priors["aux"])
  }
}

# formula tools
.generate_model_syntax             <- function(distribution, data){

  types <- attr(data[["survival"]], "type")

  model_syntax <- "model{\n"

  for(i in seq_along(types)){

    if(i == 1){
      from <- 1
      to   <- data[["survival"]][[paste0("n_", types[i])]]
    }else{
      from <- from + data[["survival"]][[paste0("n_", types[i-1])]]
      to   <- to   + data[["survival"]][[paste0("n_", types[i])]]
    }

    model_syntax <- paste0(model_syntax, .JAGS_survival_likelihood(distribution, types[i], from, to))
  }

  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}
.generate_model_formula_list       <- function(formula){

  # remove the left hand side
  formula[2] <- NULL
  formula    <- list("mu" = formula)

  return(formula)
}
.generate_model_formula_data_list  <- function(data){

  data <- list("mu" = data.frame(data[["predictors"]]))

  return(data)
}
.generate_model_formula_prior_list <- function(priors){

  priors <- list("mu" = c(priors["intercept"], priors[["terms"]]))

  return(priors)
}

# marglik tools
.marglik_function      <- function(parameters, data, distribution){

  types <- attr(data[["survival"]], "type")

  log_lik <- 0

  for(i in seq_along(types)){

    if(i == 1){
      from <- 1
      to   <- data[["survival"]][[paste0("n_", types[i])]]
    }else{
      from <- from + data[["survival"]][[paste0("n_", types[i-1])]]
      to   <- to   + data[["survival"]][[paste0("n_", types[i])]]
    }

    log_lik <- log_lik + .marglik_survival_likelihood(
      distribution = distribution,
      type         = types[type],
      x            = data[["time"]][from:to],
      mu           = parameters[["mu"]][from:to],
      aux          = parameters[["aux"]])
  }


  return(log_lik)
}
.marglik_function_null <- function(priors, data){

  types <- attr(data[["survival"]], "type")

  log_lik <- 0

  for(i in seq_along(types)){

    if(i == 1){
      from <- 1
      to   <- data[["survival"]][[paste0("n_", types[i])]]
    }else{
      from <- from + data[["survival"]][[paste0("n_", types[i-1])]]
      to   <- to   + data[["survival"]][[paste0("n_", types[i])]]
    }

    log_lik <- log_lik + .marglik_survival_likelihood(
      distribution = distribution,
      type         = types[type],
      x            = data[["time"]][from:to],
      mu           = priors[["intercept"]]$parameters[["location"]],
      aux          = priors[["aux"]]$parameters[["location"]])
  }


  return(log_lik)
}

# additional tools
.JAGS_survival_likelihood    <- function(distribution, type, from, to){
  paste0(
    "for(i in ", from,":", to, "){\n",
    "   time[i]", " ~ ", gsub("-", "_", distribution), "_", type, "(", "mu[i]", if(.has_aux(distribution)){", aux"}, ")\n",
    "}\n")
}
.marglik_survival_likelihood <- function(distribution, type, x, mu, aux){

  # TODO: add other types of censoring
  survival_likelihood <- eval(parse(text = paste0(gsub("-", "_", distribution), "_log_", switch(
    type,
    "event" = "density",
    "rcent" = "survival",
    "lcent" = NULL,
    "icent" = NULL,
    "delay" = NULL
  ))))
  args <- list(
    t   = x,
    eta = mu
  )
  if(.has_aux(distribution)){
    args <- c(args, list(aux))
  }

  log_lik <- do.call(survival_likelihood, args)

  return(log_lik)
}
.fitting_priority            <- function(models){

  terms         <- sapply(models, function(m)sum(!sapply(m[["priors"]][["terms"]], BayesTools::is.prior.point)))
  distributions <- sapply(models, function(m)switch(
    as.character(m$distribution),
    "exp-aft"     = 0,
    "weibull-aft" = 1.5,
    "lnorm-aft"   = 1,
    "llogis-aft"  = 1,
    "gamma-aft"   = 1.5
  ))

  fitting_difficulty <- terms + distributions

  return(order(fitting_difficulty, decreasing = TRUE))
}
