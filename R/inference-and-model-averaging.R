.ensemble_inference    <- function(object){

  # use only converged models with prior weights > 0 for inference about parameters
  prior_weights <- sapply(object[["models"]], function(model) model[["prior_weights"]])
  models        <- object[["models"]][.get_model_convergence(object) & prior_weights > 0]

  model_predictors      <- lapply(object[["models"]], function(model) model[["terms"]])
  model_predictors_test <- lapply(object[["models"]], function(model) model[["terms_test"]])
  model_distributions   <- sapply(object[["models"]], function(model) model[["distribution"]])

  distributions   <- object$add_info[["distributions"]]
  predictors      <- object$add_info[["predictors"]]
  predictors_test <- object$add_info[["predictors_test"]]

  # define inference options
  components      <- NULL
  parameters      <- NULL
  components_null <- list()
  parameters_null <- list()

  components_distributions      <- NULL
  components_distributions_null <- list()

  # distributions
  for(i in seq_along(distributions)){
    components_distributions                          <- c(components_distributions, distributions[i])
    components_distributions_null[[distributions[i]]] <- model_distributions != distributions[i]
  }

  # predictors
  for(i in seq_along(predictors_test)){
    components                            <- c(components, predictors_test[i])
    components_null[[predictors_test[i]]] <- sapply(model_predictors_test, function(x) !(predictors_test[i] %in% x))
  }

  for(i in seq_along(predictors)){
    parameters                                      <- c(parameters, paste0("mu_", predictors[i]))
    parameters_null[[paste0("mu_", predictors[i])]] <- sapply(model_predictors_test, function(x) !(predictors_test[i] %in% x))
  }


  ### get models inference
  inference <- BayesTools::ensemble_inference(
    model_list   = models,
    parameters   = components,
    is_null_list = components_null,
    conditional  = FALSE
  )
  # deal with the possibility of only null models models
  if(all(sapply(components_null, all))){
    inference_conditional <- NULL
  }else{
    inference_conditional <- BayesTools::ensemble_inference(
      model_list   = models,
      parameters   = components[!sapply(components_null, all)],
      is_null_list = components_null[!sapply(components_null, all)],
      conditional  = TRUE
    )
  }
  inference_distributions <- BayesTools::ensemble_inference(
    model_list   = models,
    parameters   = components_distributions[!sapply(components_distributions_null, all)],
    is_null_list = components_distributions_null[!sapply(components_distributions_null, all)],
    conditional  = FALSE
  )
  inference_distributions_conditional <- BayesTools::ensemble_inference(
    model_list   = models,
    parameters   = components_distributions[!sapply(components_distributions_null, all)],
    is_null_list = components_distributions_null[!sapply(components_distributions_null, all)],
    conditional  = TRUE
  )


  ### get model-averaged posteriors
  posteriors <- BayesTools::mix_posteriors(
    model_list   = models,
    parameters   = parameters,
    is_null_list = parameters_null,
    seed         = object$add_info[["seed"]],
    conditional  = FALSE
  )

  # deal with the possibility of only null models models
  if(all(sapply(components_null, all))){
    posteriors_conditional <- NULL
  }else{
    posteriors_conditional <- BayesTools::mix_posteriors(
      model_list   = models,
      parameters   = parameters[!sapply(parameters_null, all)],
      is_null_list = parameters_null[!sapply(parameters_null, all)],
      seed         = object$add_info[["seed"]],
      conditional  = TRUE
    )
  }

  posteriors_distributions <- .mix_posteriors_intercepts_and_aux(
    model_list          = models,
    distributions       = distributions,
    model_distributions = model_distributions,
    seed                = object$add_info[["seed"]]
  )

  # return the results
  output <- list(
    inference                           = inference,
    inference_conditional               = inference_conditional,
    inference_distributions             = inference_distributions,
    inference_distributions_conditional = inference_distributions_conditional,

    posteriors               = posteriors,
    posteriors_conditional   = posteriors_conditional,
    posteriors_distributions = posteriors_distributions
  )
  return(output)
}

.mix_posteriors_intercepts_and_aux <- function(model_list, distributions, model_distributions, seed){

  posteriors_intercepts_and_aux <- list()

  for(i in seq_along(distributions)){

    posteriors_intercepts_and_aux[paste0("(", distributions[i], ") intercept")] <- BayesTools::mix_posteriors(
      model_list   = model_list,
      parameters   = "mu_intercept",
      is_null_list = list("mu_intercept" = c(model_distributions != distributions[i])),
      seed         = seed,
      conditional  = TRUE
    )

    if(.has_aux(distributions[i])){
      posteriors_intercepts_and_aux[paste0("(", distributions[i], ") aux")] <- BayesTools::mix_posteriors(
        model_list   = model_list,
        parameters   = "aux",
        is_null_list = list("aux" = c(model_distributions != distributions[i])),
        seed         = seed,
        conditional  = TRUE
      )
    }
  }

  return(posteriors_intercepts_and_aux)
}

.compute_coeficients   <- function(RoBSA){

  return(do.call(c, unname(lapply(RoBSA$posteriors, function(posterior){
    if(inherits(posterior, "mixed_posteriors.factor")){
      out        <- apply(posterior, 2, mean)
      names(out) <- gsub("mu_", "", names(out))
    }else{
      out        <- mean(posterior)
      names(out) <- gsub("mu_", "", attr(posterior,"parameter"))
    }
    return(out)
  }))))
}
