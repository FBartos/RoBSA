#' @title Prints a fitted RoBSA object
#'
#' @param x a fitted RoBSA object.
#' @param ... additional arguments.
#'
#'
#' @return \code{print.RoBSA} invisibly returns the print statement.
#'
#' @seealso [RoBSA()]
#' @export
print.RoBSA <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimates:\n")
  print(stats::coef(x))
}


#' @title Summarize fitted RoBSA object
#'
#' @description \code{summary.RoBSA} creates a numerical
#' summary of the RoBSA object.
#'
#' @param object a fitted RoBSA object.
#' @param type whether to show the overall RoBSA results (\code{"ensemble"}),
#' an overview of the individual models (\code{"models"}), or detailed summary
#' for the individual models (\code{"individual"}).
#' @param conditional show the conditional estimates (assuming that the
#' alternative is true). Defaults to \code{FALSE}. Only available for
#' \code{type == "conditional"}.
#' @param exp whether exponents of the regression estimates should be also presented
#' @param probs quantiles of the posterior samples to be displayed.
#' Defaults to \code{c(.025, .50, .975)}
#' @param logBF show log of the BFs. Defaults to \code{FALSE}.
#' @param BF01 show BF in support of the null hypotheses. Defaults to
#' \code{FALSE}.
#' @param transform_orthonormal Whether factors with orthonormal prior
#' distributions should be transformed to differences from the grand mean. Defaults
#' to \code{TRUE}.
#' @param ... additional arguments
#' @inheritParams BayesTools::BayesTools_ensemble_tables
#'
#'
#' @return summary of a RoBSA object
#' @examples \dontrun{
#' # (execution of the example takes several minutes)
#' # example from the README (more details and explanation therein)
#' data(cancer, package = "survival")
#' priors <- calibrate_quartiles(median_t = 5, iq_range_t = 10, prior_sd = 0.5)
#' df <- data.frame(
#'   time         = veteran$time / 12,
#'   status       = veteran$status,
#'   treatment    = factor(ifelse(veteran$trt == 1, "standard", "new"), levels = c("standard", "new")),
#'   karno_scaled = veteran$karno / 100
#' )
#' RoBSA.options(check_scaling = FALSE)
#' fit <- RoBSA(
#'   Surv(time, status) ~ treatment + karno_scaled,
#'   data   = df,
#'   priors = list(
#'     treatment    = prior_factor("normal", parameters = list(mean = 0.30, sd = 0.15),
#'                                 truncation = list(0, Inf), contrast = "treatment"),
#'     karno_scaled = prior("normal", parameters = list(mean = 0, sd = 1))
#'   ),
#'   test_predictors = "treatment",
#'   prior_intercept = priors[["intercept"]],
#'   prior_aux       = priors[["aux"]],
#'   parallel = TRUE, seed = 1
#' )
#'
#' # summary can provide many details about the model
#' summary(fit)
#'
#' # note that the summary function contains additional arguments
#' # that allow to obtain a specific output, i.e, the conditional estimates
#' # (assuming that the non-null models are true) can be obtained
#' summary(fit, conditional = TRUE)
#'
#' # overview of the models and their prior and posterior probability, marginal likelihood,
#' # and inclusion Bayes factor:
#' summary(fit, type = "models")
#'
#' # and the model diagnostics overview, containing maximum R-hat and minimum ESS across parameters
#' # but see '?diagnostics' for diagnostics plots for individual model parameters
#' summary(fit, type = "diagnostics")
#'
#' # summary of individual models and their parameters can be further obtained by
#' summary(fit, type = "individual")
#'
#' }
#'
#' @note See [diagnostics()] for visual convergence checks of the individual models.
#'
#'
#' @return \code{summary.RoBSA} returns a list of tables of class 'BayesTools_table'.
#'
#' @seealso [RoBSA()], [diagnostics()], [check_RoBSA()]
#' @export
summary.RoBSA       <- function(object, type = "ensemble", conditional = FALSE,
                                exp = FALSE, parameters = FALSE, probs = c(.025, .975), logBF = FALSE, BF01 = FALSE,
                                transform_orthonormal = TRUE, short_name = FALSE, remove_spike_0 = FALSE, ...){

  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(type, "type")
  BayesTools::check_bool(exp,  "exp")
  BayesTools::check_bool(parameters,  "parameters")
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(BF01,  "BF01")
  BayesTools::check_bool(logBF, "logBF")
  BayesTools::check_bool(transform_orthonormal, "transform_orthonormal")
  BayesTools::check_bool(short_name, "short_name")
  BayesTools::check_bool(remove_spike_0, "remove_spike_0")

  # print diagnostics if all models fail to converge
  if(!any(.get_model_convergence(object))){
    if(substr(type,1,1) != "d")
      warning("All models failed to converge. Model diagnostics were printed instead.")
    type        <- "diagnostics"
  }

  if(substr(type,1,1) == "e"){

    # obtain components overview
    components_distributions <- BayesTools::ensemble_inference_table(
      inference  = object$RoBSA[["inference_distributions"]],
      parameters = names(object$RoBSA[["inference_distributions"]]),
      logBF      = logBF,
      BF01       = BF01,
      title      = "Distributions summary:"
    )

    if(length(object$RoBSA[["inference"]]) > 0){
      for(i in seq_along(object$RoBSA[["inference"]])){
        attr(object$RoBSA[["inference"]][[i]], "parameter_name") <- gsub("(mu) ", "", attr(object$RoBSA[["inference"]][[i]], "parameter_name"), fixed = TRUE)
      }
      components <- BayesTools::ensemble_inference_table(
        inference      = object$RoBSA[["inference"]],
        parameters     = names(object$RoBSA[["inference"]]),
        logBF          = logBF,
        BF01           = BF01,
        title          = "Components summary:",
      )
    }else{
      components <- NULL
    }

    # obtain estimates tables
    if(length(object$RoBSA[["posteriors"]]) > 0){
      estimates <- BayesTools::ensemble_estimates_table(
        samples               = object$RoBSA[["posteriors"]],
        parameters            = names(object$RoBSA[["posteriors"]]),
        probs                 = probs,
        title                 = "Model-averaged estimates:",
        warnings              = .collect_errors_and_warnings(object),
        transform_orthonormal = transform_orthonormal,
        formula_prefix        = FALSE
      )
      if(exp){
        estimates <- .table_add_exp(estimates)
      }
    }else{
      estimates                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
      colnames(estimates)          <- c("Mean", "Median", probs)
      class(estimates)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates))
      attr(estimates, "type")      <- rep("estimate", ncol(estimates))
      attr(estimates, "rownames")  <- TRUE
      attr(estimates, "title")     <- "Conditional estimates:"
      attr(estimates, "warnings")  <- .collect_errors_and_warnings(object)
    }


    # deal with possibly empty table in case of no alternative models
    if(is.null(object$RoBSA[["posteriors_conditional"]])){
      estimates_conditional                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
      colnames(estimates_conditional)          <- c("Mean", "Median", probs)
      class(estimates_conditional)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates_conditional))
      attr(estimates_conditional, "type")      <- rep("estimate", ncol(estimates_conditional))
      attr(estimates_conditional, "rownames")  <- TRUE
      attr(estimates_conditional, "title")     <- "Conditional estimates:"
      attr(estimates_conditional, "warnings")  <- .collect_errors_and_warnings(object)
    }else{
      estimates_conditional <- BayesTools::ensemble_estimates_table(
        samples        = object$RoBSA[["posteriors_conditional"]],
        parameters     = names(object$RoBSA[["posteriors_conditional"]]),
        probs          = probs,
        title          = "Conditional estimates:",
        warnings       = .collect_errors_and_warnings(object),
        formula_prefix = FALSE
      )
      if(exp){
        estimates_conditional <- .table_add_exp(estimates_conditional)
      }
    }

    estimates_intercept <- BayesTools::ensemble_estimates_table(
      samples    = object$RoBSA[["posteriors_intercept"]],
      parameters = names(object$RoBSA[["posteriors_intercept"]]),
      probs      = probs,
      title      = "Distribution estimates (intercept):"
    )
    if(!is.null(names(object$RoBSA[["posteriors_aux"]]))){
      estimates_aux <- BayesTools::ensemble_estimates_table(
        samples    = object$RoBSA[["posteriors_aux"]],
        parameters = names(object$RoBSA[["posteriors_aux"]]),
        probs      = probs,
        title      = "Distribution estimates (auxiliary):"
      )
    }else{
      estimates_aux                    <- data.frame(matrix(nrow = 0, ncol = length(probs) + 2))
      colnames(estimates_aux)          <- c("Mean", "Median", probs)
      class(estimates_aux)             <- c("BayesTools_table", "BayesTools_ensemble_summary", class(estimates_aux))
      attr(estimates_aux, "type")      <- rep("estimate", ncol(estimates_aux))
      attr(estimates_aux, "rownames")  <- TRUE
      attr(estimates_aux, "title")     <- "Distribution estimates (auxiliary):"
    }


    #
    # if(parameters){
    #   intercept <- cbind(intercept, "Parameter" = sapply(rownames(intercept), .intercept_name), do.call(rbind, lapply(rownames(intercept), function(distribution){
    #     do.call(.intercept_transformation(distribution), list(intercept[distribution,]))
    #   })))
    # }

    ### return results
    output <- list(
      call                     = object[["call"]],
      title                    = "Robust Bayesian survival analysis",
      components_distributions = components_distributions,
      components               = components,
      estimates                = estimates
    )

    if(conditional){
      output$estimates_conditional <- estimates_conditional
    }

    if(parameters){
      output$estimates_intercept <- estimates_intercept
      output$estimates_aux       <- estimates_aux
    }


    class(output) <- "summary.RoBSA"
    attr(output, "type") <- "ensemble"

    return(output)

  }else if(substr(type,1,1) == "m"){


    if(parameters){
      components <- list("Intercept" = "mu_intercept", "Auxiliary" = "aux")
    }else{
      components <- list()
    }

    for(i in seq_along(object$add_info[["predictors"]])){
      components[[object$add_info[["predictors"]][i]]] <- .BayesTools_parameter_name(object$add_info[["predictors"]][i])
    }

    summary <- BayesTools::ensemble_summary_table(
      models         = object[["models"]],
      parameters     = components,
      title          = "Models overview:",
      footnotes      = NULL,
      warnings       = .collect_errors_and_warnings(object),
      short_name     = short_name,
      remove_spike_0 = remove_spike_0
    )


    # add distribution column to the summary
    summary <- BayesTools::add_column(
      summary,
      column_title    = "Distribution",
      column_values   = sapply(object[["models"]], function(m) m[["distribution"]]),
      column_position = 2,
      column_type     = "string")


    output <- list(
      call       = object[["call"]],
      title      = "Robust Bayesian survival analysis",
      summary    = summary
    )

    class(output) <- "summary.RoBSA"
    attr(output, "type") <- "models"

    return(output)

  }else if(substr(type,1,1) == "d"){


    if(parameters){
      components <- list("Intercept" = "mu_intercept", "Auxiliary" = "aux")
    }else{
      components <- list()
    }

    for(i in seq_along(object$add_info[["predictors"]])){
      components[[object$add_info[["predictors"]][i]]] <- .BayesTools_parameter_name(object$add_info[["predictors"]][i])
    }

    diagnostics <- BayesTools::ensemble_diagnostics_table(
      models         = object[["models"]],
      parameters     = components,
      title          = "Diagnostics overview:",
      footnotes      = NULL,
      warnings       = .collect_errors_and_warnings(object),
      short_name     = short_name,
      remove_spike_0 = remove_spike_0
    )


    # add distribution column to the summary
    diagnostics <- BayesTools::add_column(
      diagnostics,
      column_title    = "Distribution",
      column_values   = sapply(object[["models"]], function(m) m[["distribution"]]),
      column_position = 2,
      column_type     = "string")


    output <- list(
      call        = object[["call"]],
      title       = "Robust Bayesian survival analysis",
      diagnostics = diagnostics
    )

    class(output) <- "summary.RoBSA"
    attr(output, "type") <- "diagnostics"

    return(output)

  }else if(substr(type, 1, 1) == "i"){

    output <- list(
      call       = object[["call"]],
      title      = "Robust Bayesian meta-analysis",
      models     = list()
    )

    for(i in seq_along(object[["models"]])){

      summary  <- BayesTools::model_summary_table(
        model                 = object[["models"]][[i]],
        short_name            = short_name,
        remove_spike_0        = remove_spike_0,
        formula_prefix        = FALSE
      )

      estimates <- object[["models"]][[i]][["fit_summary"]]
      attr(estimates, "warnings")  <- object[["models"]][[i]][["warnings"]]
      attr(estimates, "title")     <- "Parameter estimates:"
      if(exp){
        estimates <- .table_add_exp(estimates)
      }

      output[["models"]][[i]] <- list(
        summary   = summary,
        estimates = estimates
      )
    }

    class(output) <- "summary.RoBSA"
    attr(output, "type") <- "individual"

    return(output)

  }else{
    stop(paste0("Unknown summary type: '", type, "'."))
  }
}


.table_add_exp <- function(table){
  for(i in which(attr(table, "type") == "estimate")){
    table <- BayesTools::add_column(
      table,
      column_title    = paste0("exp(", colnames(table)[i], ")"),
      column_values   = exp(table[,i]),
      column_position = max(which(attr(table, "type") == "estimate")) + 1,
      column_type     = "estimate")
  }
  table[rownames(table) == "aux", grep("exp", colnames(table))] <- NA
  return(table)
}


#' @title Prints summary object for RoBSA method
#'
#' @param x a summary of a RoBSA object
#' @param ... additional arguments
#'
#' @return \code{print.summary.RoBSA} invisibly returns the print statement.
#'
#' @seealso [RoBSA()]
#' @export
print.summary.RoBSA <- function(x, ...){

  cat("Call:\n")
  print(x[["call"]])

  cat("\n")
  cat(x[["title"]])


  if(attr(x, "type") == "ensemble"){

    cat("\n")
    print(x[["components_distributions"]])

    if(!is.null(x[["components"]])){
      cat("\n")
      print(x[["components"]])
    }

    cat("\n")
    print(x[["estimates"]])

    if(!is.null(x[["estimates_conditional"]])){
      cat("\n")
      print(x[["estimates_conditional"]])
    }

    if(!is.null(x[["estimates_intercept"]])){
      cat("\n")
      print(x[["estimates_intercept"]])
    }

    if(!is.null(x[["estimates_aux"]])){
      cat("\n")
      print(x[["estimates_aux"]])
    }

    return(invisible())

  }else if(attr(x, "type") == "models"){

    cat("\n")
    print(x[["summary"]])

    return(invisible())

  }else if(attr(x, "type") == "diagnostics"){

    cat("\n")
    print(x[["diagnostics"]])

    return(invisible())

  }else if(attr(x, "type") == "individual"){

    for(i in seq_along(x[["models"]])){

      if(i > 1){
        cat("\n")
      }
      print(x[["models"]][[i]][["summary"]])

      cat("\n")
      print(x[["models"]][[i]][["estimates"]])
    }

    return(invisible())
  }
}

#' @title Reports whether x is a RoBSA object
#'
#' @param x an object to test
#'
#' @return is.RoBSA returns a boolean.
#' @export is.RoBSA
is.RoBSA            <- function(x){
  inherits(x, "RoBSA")
}

