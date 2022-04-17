#' @title Predict method for RoBSA objects.
#'
#' @description Predicts survival/hazard/density/mean/sd for a given
#' RoBSA object. Either predicts values for each row of a fully specified
#' \code{new_data} data.frame, or for all levels of a given \code{predictor}
#' at the mean of continuous covariate values and default factor levels or
#' covariate values specified as \code{covariates_data} data.frame.
#'
#' @param object a fitted RoBSA object
#' @param time a vector of time values at which the survival/hazard/density
#' will be predicted (for each passed data point)
#' @param new_data a data.frame containing fully specified predictors for which
#' predictions should be made
#' @param predictor an alternative input to \code{new_data} that automatically
#' generates predictions for each level of the predictor across all either across
#' levels of covariates specified by \code{covariates_data} or at the default values
#' of other predictors
#' @param covariates_data a supplementary input to \code{predictor} that specifies
#' levels of covariates for which predictions should be made
#' @param type what type of prediction should be created
#' @param summarize whether the predictions should be aggregated as mean and sd.
#' Otherwise, prediction for for posterior samples is returned.
#' @param averaged whether predictions should be combined with Bayesian model-averaging
#' or whether predictions for each individual model should be returned.
#' @param conditional whether only models assuming presence of the specified
#' \code{predictor} should be used
#' @param samples number of posterior samples to be evaluated
#' @export
predict.RoBSA <- function(object, time = NULL, new_data = NULL, predictor = NULL, covariates_data = NULL,
                          type = c("survival", "hazard", "density", "mean", "sd"),
                          summarize = TRUE, averaged = TRUE, conditional = FALSE, samples = 10000){

  predictors_all  <- attr(object$data$predictors, "variables")
  predictors_info <- attr(object$data$predictors, "variables_info")
  models          <- object[["models"]]

  BayesTools::check_real(time, "time", allow_NULL = TRUE, lower = 0, check_length = FALSE)
  BayesTools::check_char(type, "type", allow_values = c("survival", "hazard", "density", "mean", "sd"))
  if(type %in% c("survival", "hazard", "density") && is.null(time))
    stop("Time variable 'time' must be provided.")
  BayesTools::check_char(predictor, "predictor", allow_NULL = TRUE, allow_values = predictors_all)
  BayesTools::check_bool(summarize, "summarize")
  BayesTools::check_bool(averaged, "averaged")
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_int(samples, "samples", lower = 0)


  # check that the new data are correctly specified
  if(!is.null(new_data) && (is.null(predictor) && is.null(covariates_data)) ){

    if(!is.data.frame(new_data))
      stop("'new_data' must be a data.frame")
    if(!all(all_predictors %in% colnames(new_data)))
      stop("All predictors must be provided.")

  }else if(is.null(new_data)){

    if(!is.null(covariates_data) && !is.data.frame(covariates_data))
      stop("'covariates_data' must be a data.frame")

    missing_predictors <- predictors_all[!predictors_all %in% c(colnames(covariates_data), predictor)]

    # add the predictor if specified
    if(!is.null(predictor)){

      predictor_data  <- data.frame(
        if(predictors_info[[predictor]][["type"]] == "factor") predictors_info[[predictor]][["levels"]]
        else if(predictors_info[[predictor]][["type"]] == "continuous") predictors_info[[predictor]][["mean"]]
      )

      if(!is.null(covariates_data)){
        covariates_data <- do.call(rbind, lapply(1:nrow(predictor_data), function(i) covariates_data))
        predictor_data  <- do.call(rbind, lapply(1:nrow(covariates_data), function(i) predictor_data))
        new_data        <- cbind(covariates_data, predictor_data)
      }else{
        new_data        <- predictor_data
      }

      colnames(new_data)[ncol(new_data)] <- predictor

    }

    # construct the missing predictors (i.e., the first level of factors and mean of covariates)
    for(i in seq_along(missing_predictors)){
      new_data <- cbind(
        new_data,
        if(predictors_info[[missing_predictors[i]]][["type"]] == "factor") predictors_info[[missing_predictors[i]]][["default"]]
        else if(predictors_info[[missing_predictors[i]]][["type"]] == "continuous") predictors_info[[missing_predictors[i]]][["mean"]])
      colnames(new_data)[ncol(new_data)] <- missing_predictors[i]
    }

  }else{
    stop("Either no data or only the 'new_data' or 'parameter' (and 'covariates_data') need to be specified.")
  }


  # rescale the new data if the original input was re-scaled
  if(object[["add_info"]][["rescale_data"]]){
    for(i in seq_along(predictors_all)){
      if(predictors_info[[predictors_all[i]]][["type"]] == "continuous"){
        new_data[,predictors_all[i]] <- .pred_scale(new_data[,predictors_all[i]], predictors_info[[predictors_all[i]]])
      }
    }
  }


  # select conditional models (if the predictor is to be tested)
  if(conditional){
    if(!predictor %in% object$add_info[["predictors_test"]])
      stop("Conditional models cannot be selected since the given predictor was not tested (i.e., it was included in all models). ")

    is_null_models <- attr(object$RoBSA$inference_conditional[[.BayesTools_parameter_name(predictor)]], "is_null")
    models         <- models[!is_null_models]
  }


  # pre-specify number of samples for model-averaging
  if(averaged){
    # get re-standardized model weights in case a condition models were selected
    model_weights <- sapply(models, function(model) model[["inference"]][["post_prob"]])
    model_weights <- model_weights / sum(model_weights)

    samples <- round(model_weights * samples)
  }else{
    samples <- rep(samples, length(models))
  }


  # obtain evaluated posterior distributions from each model
  model_parameters <- lapply(seq_along(models), function(i){

    # the samples
    posteriors <- list(
        mu = BayesTools::JAGS_evaluate_formula(
          fit        = models[[i]][["fit"]],
          formula    = object[["formula"]],
          parameter  = "mu",
          data       = new_data,
          prior_list = attr(models[[i]][["fit"]], "prior_list")
    ))
    if(.has_aux(models[[i]][["distribution"]])){
      posteriors$aux = .extract_aux_samples(models[[i]][["fit"]])
    }

    # subset them to have equal amount across models
    indx <- sample(ncol(posteriors[["mu"]]), size = samples[i], replace = TRUE)
    posteriors[["mu"]] <- posteriors[["mu"]][,indx]
    if(.has_aux(models[[i]][["distribution"]])){
      posteriors[["aux"]] <- posteriors[["aux"]][indx]
    }

    attr(posteriors, "model")        <- models[[i]][["inference"]][["m_number"]]
    attr(posteriors, "distribution") <- models[[i]][["distribution"]]

    return(posteriors)
  })


  # deal with different types of predictions separately
  if(type %in% c("mean", "sd")){

    model_predictions <- lapply(model_parameters, function(parameters){

      # get the corresponding function
      prediction_function <- eval(parse(text = paste0(gsub("-", "_", attr(parameters, "distribution")), "_", type)))

      # make predictions for each data point and merge into a data.frame
      out <- data.frame(do.call(rbind, lapply(1:nrow(new_data), function(i){
        args <- list(eta = parameters[["mu"]][i,])
        if(.has_aux(attr(parameters, "distribution"))){
          args <- c(args, list(parameters[["aux"]]))
        }
        return(do.call(prediction_function, args))
      })))

      attr(out, "data")         <- new_data
      attr(out, "model")        <- attr(parameters, "model")
      attr(out, "distribution") <- attr(parameters, "distribution")
      attr(out, "outcome")      <- type
      return(out)
    })


    # return the properly aggregated results
    if(!averaged && !summarize){

      attr(model_predictions, "data")    <- new_data
      attr(model_predictions, "outcome") <- type
      return(model_predictions)

    }else if(!averaged && summarize){

      model_predictions <- lapply(model_predictions, function(predictions){

        if(anyNA(predictions))
          warning("Some of the model returned undefined predictions. The returned summary ommits all NA samples and might be biased.", immediate. = TRUE)

        out <- data.frame(
          mean   = apply(predictions, 1, mean, na.rm = TRUE),
          sd     = apply(predictions, 1, sd,   na.rm = TRUE),
          lCI    = apply(predictions, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
          median = apply(predictions, 1, stats::quantile, probs = 0.500, na.rm = TRUE),
          uCI    = apply(predictions, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
        )

        attr(out, "data")         <- new_data
        attr(out, "model")        <- attr(predictions, "model")
        attr(out, "distribution") <- attr(predictions, "distribution")
        attr(out, "outcome")      <- type
        return(out)
      })

      attr(model_predictions, "data")    <- new_data
      attr(model_predictions, "outcome") <- type
      return(model_predictions)

    }else{

      # average the predictions
      data_predictions <- do.call(cbind, lapply(seq_along(model_predictions), function(i){
        if(samples[i] > 0){
         return(model_predictions[[i]][,1:samples[i]])
        }else{
          return(matrix(nrow = nrow(model_predictions[[i]]), ncol = 0))
        }
      }))
      colnames(data_predictions) <- NULL

      if(summarize){

        if(anyNA(data_predictions))
          warning("Some of the model returned undefined predictions. The returned summary ommits all NA samples and might be biased.", immediate. = TRUE)

        data_predictions <- data.frame(
          mean   = apply(data_predictions, 1, mean, na.rm = TRUE),
          sd     = apply(data_predictions, 1, sd,   na.rm = TRUE),
          lCI    = apply(data_predictions, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
          median = apply(data_predictions, 1, stats::quantile, probs = 0.500, na.rm = TRUE),
          uCI    = apply(data_predictions, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
        )
      }

      attr(data_predictions, "data")    <- new_data
      attr(data_predictions, "outcome") <- type
      return(data_predictions)
    }
  }

  if(type %in% c("survival", "hazard", "density")){

    model_predictions <- lapply(model_parameters, function(parameters){

      # get the corresponding function
      prediction_function <- eval(parse(text = paste0(gsub("-", "_", attr(parameters, "distribution")), "_", type)))

      # at each data point, make predictions for all times for each posterior sample
      out <- lapply(1:nrow(new_data), function(i){
        out <- do.call(rbind, lapply(seq_along(time), function(t){
          args <- list(
            t   = time[t],
            eta = parameters[["mu"]][i,]
          )
          if(.has_aux(attr(parameters, "distribution"))){
            args <- c(args, list(parameters[["aux"]]))
          }
          return(do.call(prediction_function, args))
        }))

        attr(out, "time") <- time
        return(out)
      })

      attr(out, "data")         <- new_data
      attr(out, "model")        <- attr(parameters, "model")
      attr(out, "distribution") <- attr(parameters, "distribution")
      attr(out, "outcome")      <- type
      return(out)
    })


    # return the properly aggregated results
    if(!averaged && !summarize){

      attr(model_predictions, "data")    <- new_data
      attr(model_predictions, "time")    <- time
      attr(model_predictions, "outcome") <- type
      return(model_predictions)

    }else if(!averaged && summarize){

      model_predictions <- lapply(model_predictions, function(predictions){

        if(anyNA(predictions))
          warning("Some of the model returned undefined predictions. The returned summary ommits all NA samples and might be biased.", immediate. = TRUE)

        out <- lapply(seq_along(predictions), function(i){
          out <- data.frame(
            time   = time,
            mean   = apply(predictions[[i]], 1, mean, na.rm = TRUE),
            sd     = apply(predictions[[i]], 1, sd,   na.rm = TRUE),
            lCI    = apply(predictions[[i]], 1, stats::quantile, probs = 0.025, na.rm = TRUE),
            median = apply(predictions[[i]], 1, stats::quantile, probs = 0.500, na.rm = TRUE),
            uCI    = apply(predictions[[i]], 1, stats::quantile, probs = 0.975, na.rm = TRUE)
          )
        })

        attr(out, "data")         <- attr(predictions,"data")
        attr(out, "model")        <- attr(predictions,"model")
        attr(out, "distribution") <- attr(predictions,"distribution")
        return(out)
      })

      attr(model_predictions, "data")    <- new_data
      attr(model_predictions, "time")    <- time
      attr(model_predictions, "outcome") <- type
      return(model_predictions)

    }else{

      # average the predictions
      data_predictions <- lapply(1:nrow(new_data), function(j){

        out <- do.call(cbind, lapply(seq_along(model_predictions), function(i){
          if(samples[i] > 0){
            return(model_predictions[[i]][[j]][,1:samples[i]])
          }else{
            return(matrix(nrow = nrow(model_predictions[[i]][[j]]), ncol = 0))
          }
        }))

        attr(out, "time") <- time
        return(out)
      })


      if(summarize){
        data_predictions <- lapply(data_predictions, function(predictions){

          if(anyNA(predictions))
            warning("Some of the model returned undefined predictions. The returned summary ommits all NA samples and might be biased.", immediate. = TRUE)

          return(data.frame(
            time   = time,
            mean   = apply(predictions, 1, mean, na.rm = TRUE),
            sd     = apply(predictions, 1, sd,   na.rm = TRUE),
            lCI    = apply(predictions, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
            median = apply(predictions, 1, stats::quantile, probs = 0.500, na.rm = TRUE),
            uCI    = apply(predictions, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
          ))
        })
      }

      attr(data_predictions, "data")    <- new_data
      attr(data_predictions, "time")    <- time
      attr(data_predictions, "outcome") <- type
      return(data_predictions)
    }
  }
}
