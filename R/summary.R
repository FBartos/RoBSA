#' @title Prints a fitted RoBSA object
#'
#' @param x a fitted RoBSA object.
#' @param ... additional arguments.
#' @export  print.RoBSA
#' @rawNamespace S3method(print, RoBSA)
#' @seealso [RoBSA()]
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
#' @param diagnostics show the maximum R-hat and minimum ESS for the main
#' parameters in each of the models. Only available for \code{type = "ensemble"}.
#' @param include_theta whether the estimated random effects should be included
#' either in the summaries.
#' @param probs quantiles of the posterior samples to be displayed.
#' Defaults to \code{c(.025, .50, .975)}
#' @param logBF show log of the BFs. Defaults to \code{FALSE}.
#' @param BF01 show BF in support of the null hypotheses. Defaults to
#' \code{FALSE}.
#' @param digits_estimates a number of decimals for rounding the estimates.
#' Defaults to \code{3}.
#' @param digits_BF a number of decimals for rounding the BFs. Defaults to \code{3}.
#' @param ... additional arguments
#'
#' @return summary of a RoBSA object
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBSA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
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
#' summary(fit, type = "models", diagnostics = TRUE)
#'
#' # summary of individual models and their parameters can be further obtained by
#' summary(fit, type = "individual")
#'
#' }
#' @note See [diagnostics()] for visual convergence checks of the individual models.
#' @method summary RoBSA
#' @export summary.RoBSA
#' @rawNamespace S3method(summary, RoBSA)
#' @seealso [RoBSA()] [diagnostics()]
summary.RoBSA       <- function(object, type = "ensemble", conditional = FALSE,
                                exp = FALSE, parameters = FALSE, probs = c(.025, .975), logBF = FALSE, BF01 = FALSE,
                                short_name = FALSE, remove_spike_0 = FALSE, ...){

  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(type, "type")
  BayesTools::check_bool(exp,  "exp")
  BayesTools::check_bool(parameters,  "parameters")
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  BayesTools::check_bool(BF01,  "BF01")
  BayesTools::check_bool(logBF, "logBF")
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

    components <- BayesTools::ensemble_inference_table(
      inference  = object$RoBSA[["inference"]],
      parameters = names(object$RoBSA[["inference"]]),
      logBF      = logBF,
      BF01       = BF01,
      title      = "Components summary:"
    )

    # obtain estimates tables
    estimates <- BayesTools::ensemble_estimates_table(
      samples    = object$RoBSA[["posteriors"]],
      parameters = names(object$RoBSA[["posteriors"]]),
      probs      = probs,
      title      = "Model-averaged estimates:",
      warnings   = .collect_errors_and_warnings(object)
    )

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
        samples    = object$RoBSA[["posteriors_conditional"]],
        parameters = names(object$RoBSA[["posteriors_conditional"]]),
        probs      = probs,
        title      = "Conditional estimates:",
        warnings   = .collect_errors_and_warnings(object)
      )
    }

    estimates_intercept <- BayesTools::ensemble_estimates_table(
      samples    = object$RoBSA[["posteriors_intercept"]],
      parameters = names(object$RoBSA[["posteriors_intercept"]]),
      probs      = probs,
      title      = "Distribution estimates (intercept):"
    )
    estimates_aux <- BayesTools::ensemble_estimates_table(
      samples    = object$RoBSA[["posteriors_aux"]],
      parameters = names(object$RoBSA[["posteriors_aux"]]),
      probs      = probs,
      title      = "Distribution estimates (auxiliary):"
    )


    ### return results
    output <- list(
      call                     = object[["call"]],
      title                    = "Robust Bayesian meta-analysis",
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
      components <- list("Intercept" = "mu_intercept", "Axillary" = "aux")
    }else{
      components <- list()
    }

    for(i in seq_along(object$add_info[["predictors"]])){
      components[[object$add_info[["predictors"]][i]]] <- paste0("mu_", object$add_info[["predictors"]][i])
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

    .BayesTools_table_add_distributions(summary, sapply(object[["models"]], function(m) m[["distribution"]]))
    # add distribution



    output <- list(
      call       = object[["call"]],
      title      = "Robust Bayesian meta-analysis",
      summary    = summary
    )

    class(output) <- "summary.RoBSA"
    attr(output, "type") <- "models"

    return(output)

  }else if(substr(type,1,1) == "d"){

    components <- names(object$RoBSA[["inference"]])[names(object$RoBSA[["inference"]]) %in% c("Effect", "Heterogeneity", "Bias")]
    parameters <- list()
    if(any(components == "Effect")){
      parameters[["Effect"]] <- "mu"
    }
    if(any(components == "Heterogeneity")){
      parameters[["Heterogeneity"]] <- "tau"
      if(!attr(object$data, "all_independent")){
        parameters[["Var. allocation"]] <- "rho"
      }
    }
    if(any(components == "Bias")){
      parameters[["Bias"]] <- c("PET", "PEESE", "omega")
    }

    diagnostics <- BayesTools::ensemble_diagnostics_table(
      models         = object[["models"]],
      parameters     = parameters,
      title          = "Diagnostics overview:",
      footnotes      = NULL,
      warnings       = .collect_errors_and_warnings(object),
      short_name     = short_name,
      remove_spike_0 = remove_spike_0
    )

    output <- list(
      call        = object[["call"]],
      title       = "Robust Bayesian meta-analysis",
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
        model          = object[["models"]][[i]],
        short_name     = short_name,
        remove_spike_0 = remove_spike_0
      )
      if(output_scale == "y"){
        estimates <- object[["models"]][[i]][["fit_summary"]]
        attr(estimates, "warnings")  <- object[["models"]][[i]][["warnings"]]
        attr(estimates, "title")     <- "Parameter estimates:"
      }else{
        estimates <- object[["models"]][[i]][["fit_summaries"]][[output_scale]]
        attr(estimates, "footnotes") <- .scale_note(object[["models"]][[i]][["prior_scale"]], output_scale)
        attr(estimates, "warnings")  <- object[["models"]][[i]][["warnings"]]
        attr(estimates, "title")     <- "Parameter estimates:"
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





  if(substr(type,1,1) == "e"){

    ### model estimates
    predictors <- list()
    intercept  <- list()
    aux        <- list()

    # compute quantiles
    if(!is.null(probs)){
      if(length(probs) != 0){

        if(!is.numeric(probs) | !is.vector(probs))stop("The passed probabilities 'probs' must be a numeric vector.")
        if(!(all(probs > 0) & all(probs < 1)))stop("The passed probabilities 'probs' must be higher than 0 and lower than 1.")

        # quantiles
        for(type in c("averaged", "conditional")){
          if(is.null(attr(object$data, "predictors"))){
            predictors[[type]] <- data.frame(matrix(ncol = 2, nrow = 0))
            colnames(predictors[[type]]) <- probs
          }else{
            predictors[[type]] <- data.frame(do.call(rbind, lapply(object$RoBSA$samples[["predictors"]][[type]], stats::quantile, probs = probs)))
            colnames(predictors[[type]]) <- probs
          }

        }
        intercept <- data.frame(do.call(rbind, lapply(object$RoBSA$samples[["intercept"]], stats::quantile, probs = probs)))
        aux       <- data.frame(do.call(rbind, lapply(object$RoBSA$samples[["aux"]],       stats::quantile, probs = probs)))
        colnames(intercept) <- probs
        colnames(aux)      <- probs
      }
    }

    # pointe stimates
    for(type in c("averaged", "conditional")){
      if(is.null(attr(object$data, "predictors"))){
        predictors[[type]] <- cbind("Mean" = numeric(), "Median" = numeric(), predictors[[type]])
      }else{
        predictors[[type]] <- cbind(do.call(rbind, lapply(object$RoBSA$samples[["predictors"]][[type]], function(x)c("Mean" = mean(x), "Median" = stats::median(x)))), predictors[[type]])
      }

    }
    intercept <- cbind(do.call(rbind, lapply(object$RoBSA$samples[["intercept"]], function(x)c("Mean" = mean(x), "Median" = stats::median(x)))), intercept)
    aux       <- cbind(do.call(rbind, lapply(object$RoBSA$samples[["aux"]],       function(x)c("Mean" = mean(x), "Median" = stats::median(x)))), aux)


    ### adding transformed results
    if(exp){
      for(type in c("averaged", "conditional")){
        predictors[[type]] <- cbind(predictors[[type]], exp(predictors[[type]]))
        colnames(predictors[[type]])[(ncol(predictors[[type]])/2+1):ncol(predictors[[type]])] <-
          paste0("exp(",colnames(predictors[[type]])[(ncol(predictors[[type]])/2+1):ncol(predictors[[type]])],")")
      }
    }

    if(parameters){
      intercept <- cbind(intercept, "Parameter" = sapply(rownames(intercept), .intercept_name), do.call(rbind, lapply(rownames(intercept), function(distribution){
        do.call(.intercept_transformation(distribution), list(intercept[distribution,]))
      })))
    }


    ### adding tests to the predictors
    if(is.null(attr(object$data, "predictors"))){
      predictors[["averaged"]]  <- cbind(predictors[["averaged"]], numeric(), numeric(), numeric())
      colnames(predictors[["averaged"]])[(ncol(predictors[["averaged"]])-2):ncol(predictors[["averaged"]])]  <- c(
        "Prior prob.",
        "Post. prob.",
        paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")")
      )
    }else{
      predictors_prior_prob     <- unlist(object$RoBSA$prior_prob$predictors)
      predictors_posterior_prob <- unlist(object$RoBSA$posterior_prob$predictors)
      predictors_BF             <- unlist(object$RoBSA$BF$predictors)
      predictors_BF             <- .BF_format(predictors_BF, BF01, logBF)
      predictors[["averaged"]]  <- cbind(predictors[["averaged"]], predictors_prior_prob, predictors_posterior_prob, predictors_BF)
      colnames(predictors[["averaged"]])[(ncol(predictors[["averaged"]])-2):ncol(predictors[["averaged"]])]  <- c(
        "Prior prob.",
        "Post. prob.",
        paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")")
      )
    }



    ### fixing naming
    rownames(aux) <- paste0(sapply(rownames(aux), .aux_name), " (", rownames(aux), ")")


    ### model types overview
    distributons                <- as.character(sapply(object$models, function(m)m$distribution))
    distributons_prior_prob     <- unlist(object$RoBSA$prior_prob$distributions)
    distributons_posterior_prob <- unlist(object$RoBSA$posterior_prob$distributions)
    distributons_BF             <- unlist(object$RoBSA$BF$distributions)
    distributons_BF             <- .BF_format(distributons_BF, BF01, logBF)
    distributons                <- distributons[distributons %in% names(distributons_BF)]

    overview_tab <- data.frame(
      n_models   = as.vector(table(distributons)[unique(distributons)]),
      prior_prob = distributons_prior_prob[unique(distributons)],
      post_prob  = distributons_posterior_prob[unique(distributons)],
      BF         = distributons_BF[unique(distributons)]
    )
    colnames(overview_tab) <- c(
      "Models",
      "Prior prob.",
      "Post. prob.",
      paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")")
    )

    if(!conditional){
      predictors[["conditional"]] <- NULL
    }
    if(!parameters){
      aux <- NULL
    }

    ### return results
    res <- list(
      call       = object$call,
      overview   = overview_tab,
      predictors = predictors,
      intercept  = intercept,
      aux        = aux,
      add_info   = list(
        n_models         = length(distributons),
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        type             = "ensemble",
        failed           = sum(!object$add_info$converged)
      )
    )


  }else if(substr(type,1,1) == "m"){

    predictors        <- attr(object$data, "predictors")
    distibutions      <- sapply(1:length(object$models), function(i)object$models[[i]]$distribution)

    priors_intercept  <- sapply(1:length(object$models), function(i)print(object$models[[i]]$priors[["intercept"]], silent = TRUE))
    priors_aux        <- sapply(1:length(object$models), function(i){
      if(is.null(object$models[[i]]$priors[["aux"]])){
        return("")
      }else{
        return(print(object$models[[i]]$priors[["aux"]], silent = TRUE))
      }
    })
    priors_predictors <- lapply(1:length(object$models), function(i){
      priors <- NULL
      for(j in seq_along(predictors)){
        priors <- c(priors, print(object$models[[i]]$priors[["predictors"]][[predictors[j]]], silent = TRUE))
      }
      return(priors)
    })
    if(!is.null(predictors)){
      priors_predictors <- do.call(rbind, priors_predictors)
    }
    prior_odds     <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_odds)
    prior_prob     <- prior_odds / sum(prior_odds)
    marg_lik       <- sapply(1:length(object$models), function(i)object$models[[i]]$marg_lik$logml)
    posterior_prob <- bridgesampling::post_prob(marg_lik, prior_prob = prior_prob)
    BF             <- sapply(1:length(object$models), function(i){
      temp_mm <- rep(FALSE, length(object$models))
      temp_mm[i] <- TRUE
      .inclusion_BF(prior_prob, posterior_prob, temp_mm)
    })

    BF[is.nan(BF)] <- NA
    if(BF01){
      BF <- 1/BF
    }else{
      BF <- BF
    }
    if(logBF){
      BF <- log(BF)
    }

    overview_tab <- data.frame("Distribution" = distibutions)
    for(i in seq_along(predictors)){
      overview_tab <- cbind(overview_tab, priors_predictors[,i])
      colnames(overview_tab)[ncol(overview_tab)] <- paste0("Prior ", predictors[i])
    }

    overview_tab2 <- data.frame(
      priors_intercept,
      priors_aux,
      prior_prob,
      posterior_prob,
      marg_lik,
      BF,
      stringsAsFactors = FALSE
    )
    rownames(overview_tab2) <- NULL
    colnames(overview_tab2) <- c("Prior intercept", "Prior aux", "Prior prob.", "Post. prob.", "log(MargLik)",
                                paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")"))

    overview_tab <- cbind(overview_tab, overview_tab2)
    # add the summary model diagnostics
    if(diagnostics){

      diagnostics_tab <- overview_tab[,c(1, grep("Prior", colnames(overview_tab)))]
      diagnostics_tab <- diagnostics_tab[,-ncol(diagnostics_tab)]

      # extract max(R-hat) & min(ESS)
      diag_sum <- sapply(1:length(object$models), function(i){

        temp_x <- object$models[[i]]$fit

        if(length(temp_x) == 0 | any(class(object$models[[i]]$fit) %in% c("simpleError","error"))){

          return(c(NA, NA, NA))

        }else{

          s.x <- object$models[[i]]$fit_summary

          if(length(dim(s.x)) == 2){
            return(c(
              MCerr = max(s.x[, 7]),
              Rhat  = max(s.x[,11]),
              ESS   = min(s.x[, 9])
            ))
          }else{
            return(c(
              MCerr = max(s.x[ 7]),
              Rhat  = max(s.x[11]),
              ESS   = min(s.x[ 9])
            ))
          }

        }
       })

      diagnostics_tab$"max(MCMC error)"  <- diag_sum[1,]
      diagnostics_tab$"min(ESS)"         <- diag_sum[3,]
      diagnostics_tab$"max(Rhat)"        <- diag_sum[2,]
    }


    res <- list(
      call     = object$call,
      overview = overview_tab,
      add_info = list(
        n_models         = length(object$models),
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        predictors       = predictors,
        distibutions     = distibutions,
        type             = "models"
      )
    )

    if(diagnostics){
      res$diagnostics <- diagnostics_tab
    }


  }else if(substr(type, 1, 1) == "i"){

    overview_tabs <- list()

    prior_odds     <- sapply(1:length(object$models), function(i)object$models[[i]]$prior_odds)
    prior_prob     <- prior_odds / sum(prior_odds)
    posterior_prob <- object$RoBSA$posterior_prob$all
    marg_lik       <- sapply(1:length(object$models), function(i)object$models[[i]]$marg_lik$logml)
    BF             <- sapply(1:length(object$models), function(i){
      temp_mm <- rep(FALSE, length(object$models))
      temp_mm[i] <- TRUE
      .inclusion_BF(prior_prob, posterior_prob, temp_mm)
    })
    BF             <- .BF_format(BF, BF01, logBF)


    for(i in 1:length(object$models)){

      if(length(object$models[[i]]$fit) == 0 |  any(class(object$models[[i]]$fit) %in% c("simpleError","error"))){

        s.x <- NULL

      }else{

        s.x <- object$models[[i]]$fit_summary[,c(4, 5, 1, 2, 3, 7, 8, 9, 11), drop = FALSE]

        colnames(s.x) <- c("Mean", "SD", ".025", "Median", ".975", "MCMC error", "Error % of SD", "ESS", "Rhat")

      }

      s.x_predictors  <- s.x[grepl("beta",  rownames(s.x)),,drop = FALSE]
      s.x_int_and_aux <- s.x[!grepl("beta", rownames(s.x)),,drop = FALSE]

      ### adding transformed results
      if(!is.null(s.x_predictors)){
        if(exp){
          s.x_predictors <- cbind(s.x_predictors[,1:5, drop = FALSE], exp(s.x_predictors[,c(1,3,4,5), drop = FALSE]), s.x_predictors[,6:9, drop = FALSE])
          colnames(s.x_predictors)[6:9] <- paste0("exp(",colnames(s.x_predictors)[6:9],")")
        }
      }

      if(!is.null(s.x_predictors)){
        if(parameters){
          s.x_int_and_aux[rownames(s.x_int_and_aux) == "intercept",1:5] <- do.call(.intercept_transformation(as.character(object$models[[i]]$distribution)), list(s.x_int_and_aux[rownames(s.x_int_and_aux) == "intercept",1:5]))
          s.x_int_and_aux[rownames(s.x_int_and_aux) == "intercept",2] <- NA
          rownames(s.x_int_and_aux)[rownames(s.x_int_and_aux) == "intercept"] <- .intercept_name(as.character(object$models[[i]]$distribution))
        }
        if(.has_aux(as.character(object$models[[i]]$distribution))){
          rownames(s.x_int_and_aux)[rownames(s.x_int_and_aux) == "aux"] <- .aux_name(as.character(object$models[[i]]$distribution))
        }
      }

      overview_tabs[[i]] <- list(
        priors          = object$models[[i]]$priors,
        tab_predictors  = s.x_predictors,
        tab_int_and_aux = s.x_int_and_aux,
        prior_prob      = prior_prob[i],
        marg_lik        = marg_lik[i],
        posterior_prob  = posterior_prob[i],
        BF              = BF[i],
        distribution    = as.character(object$models[[i]]$distribution),
        add_info        = list(
          BF_type      = paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")")
        )
      )
    }


    res <- list(
      call     = object$call,
      overview = overview_tabs,
      add_info = list(
        n_models         = length(object$models),
        digits_estimates = digits_estimates,
        digits_BF        = digits_BF,
        exp              = exp,
        parameters       = parameters,
        BF_type          = paste0("Incl. ", if(logBF)"log(",if(BF01)"1/","BF",if(logBF)")"),
        type             = "individual"
      )

    )

  }

  class(res) <- "summary.RoBSA"
  return(res)
}


.BayesTools_table_add_distributions <- function(summary, distributions){

  new_summary <- cbind("Model" = summary[,1], "Distribution" = sapply(object[["models"]], function(m) m[["distribution"]]), summary[,-1])

  attr(new_summary, "type") <- c(attr(summary, "type")[1], "string", attr(summary, "type")[-1])

  for(a in names(attributes(summary))[!names(attributes(summary)) %in% c("names", "row.names", "class", "type")]){
    attr(new_summary, a) <- attr(summary, a)
  }

  return(new_summary)
}


#' @title Prints summary object for RoBSA method
#'
#' @param x a summary of a RoBSA object
#' @param ... additional arguments
#' @method print.summary RoBSA
#' @export print.summary.RoBSA
#' @rawNamespace S3method(print, summary.RoBSA)
#' @seealso [RoBSA()]
print.summary.RoBSA <- function(x, ...){

  # format the output before printing
  if(x$add_info$type == "ensemble"){

    overview <- x$overview
    overview$Models <- paste0(overview$Models,"/",  x$add_info$n_models - x$add_info$failed)
    overview[,2:3]  <- format(round(overview[,2:3], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
    overview[,4]    <- format(round(overview[,4],   x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)

    # round the results (a loop is required to deal with NAs)
    averaged <- x$predictors$averaged
    for(i in 1:ncol(averaged)){
      if(i == ncol(averaged)){
        averaged[,i] <- format(round(averaged[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }else{
        averaged[,i] <- format(round(averaged[,i], x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)
      }
    }
    if(!is.null(x$predictors$conditional)){
      conditional <- x$predictors$conditional
      for(i in 1:ncol(conditional)){
        conditional[,i] <- format(round(conditional[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }
    }
    if(!is.null(x$intercept)){
      intercept <- x$intercept
      for(i in 1:ncol(intercept)){
        if(colnames(intercept)[i] != "Parameter"){
          intercept[,i] <- format(round(intercept[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
        }
      }
    }
    if(!is.null(x$aux)){
      aux <- x$aux
      for(i in 1:ncol(aux)){
        aux[,i] <- format(round(aux[,i], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }
    }

  }else if(x$add_info$type == "models"){

    overview <- x$overview
    for(cn in c("Prior prob.",  "Post. prob.", "log(MargLik)")){
      overview[,cn]  <- format(round(overview[,cn], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
    }
    overview[,ncol(overview)]    <- format(round(overview[,ncol(overview)],   x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)

    if(!is.null(x$diagnostics)){
      diagnostics     <- x$diagnostics
      for(cn in c("max(MCMC error)",  "max(Rhat)")){
        diagnostics[,cn]  <- format(round(diagnostics[,cn], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      }
      diagnostics[,"min(ESS)"]    <- format(round(diagnostics[,"min(ESS)"],   x$add_info$digits_BF),        nsmall = x$add_info$digits_BF)
    }
  }else if(x$add_info$type == "individual"){

    overview <- x$overview

    for(i in 1:length(overview)){

      temp_main_info_names <- c("Model:", "Distribution", "Prior prob.:", "log(MargLik):", "Post. prob.:", paste0(x$add_info$BF_type,":"))
      temp_main_info       <- c(i, overview[[i]]$distribution,
                                format(round(overview[[i]]$prior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates),
                                format(round(overview[[i]]$marg_lik, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates),
                                format(round(overview[[i]]$posterior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates),
                                format(round(overview[[i]]$BF, x$add_info$digits_BF), nsmall = x$add_info$digits_BF))


      temp_main_priors_names <- "Parameter"
      temp_main_priors       <- "Prior Distributions"
      if(!is.null(overview[[i]]$priors$predictors)){
        for(j in seq_along(overview[[i]]$priors$predictors)){
          temp_main_priors_names <- c(temp_main_priors_names, names(overview[[i]]$priors$predictors)[j])
          temp_main_priors       <- c(temp_main_priors, print(overview[[i]]$priors$predictors[[j]], silent = TRUE))
        }
      }
      temp_main_priors_names <- c(temp_main_priors_names, "intercept")
      temp_main_priors       <- c(temp_main_priors, print(overview[[i]]$priors$intercept, silent = TRUE))
      if(.has_aux(overview[[i]]$distribution)){
        temp_main_priors_names <- c(temp_main_priors_names, .aux_name(overview[[i]]$distribution))
        temp_main_priors       <- c(temp_main_priors, print(overview[[i]]$priors$aux, silent = TRUE))
      }

      if(length(temp_main_info) > length(temp_main_priors)){
        temp_main_priors       <- c(temp_main_priors,       rep("", length(temp_main_info) - length(temp_main_priors)))
        temp_main_priors_names <- c(temp_main_priors_names, rep("", length(temp_main_info) - length(temp_main_priors_names)))
      }else if(length(temp_main_priors) > length(temp_main_info)){
        temp_main_info       <- c(temp_main_priors,       rep("", length(temp_main_priors) - length(temp_main_info)))
        temp_main_info_names <- c(temp_main_priors_names, rep("", length(temp_main_priors) - length(temp_main_info_names)))
      }

      temp_main <- data.frame(cbind(
        temp_main_info_names, temp_main_info,
        rep("", length(temp_main_info)), rep("", length(temp_main_info)),
        temp_main_priors_names, temp_main_priors))
      colnames(temp_main) <- c("", " ", "  ", "   ", "    ", "     ")
      overview[[i]]$tab_main <- temp_main

      ind_pred  <- if(x$add_info$exp) c(1:10,13)  else c(1:6,9)
      ind_pred2 <- if(x$add_info$exp) 11          else 7
      ind_pred3 <- if(x$add_info$exp) 12          else 8

      if(!is.null(overview[[i]]$tab_predictors)){
        temp_tab <- overview[[i]]$tab_predictors
        overview[[i]]$tab_predictors[,ind_pred]  <- format(round(temp_tab[,ind_pred], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
        overview[[i]]$tab_predictors[,ind_pred2] <- format(round(temp_tab[,ind_pred2], 1), nsmall = 1)
        overview[[i]]$tab_predictors[,ind_pred3] <- round(temp_tab[,ind_pred3])
      }

      if(!is.null(overview[[i]]$tab_int_and_aux)){
        temp_tab <- overview[[i]]$tab_int_and_aux
        overview[[i]]$tab_int_and_aux[,c(1:6,9)]  <- format(round(temp_tab[,c(1:6,9)], x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
        overview[[i]]$tab_int_and_aux[,7]         <- format(round(temp_tab[,7], 1), nsmall = 1)
        overview[[i]]$tab_int_and_aux[,8]         <- round(temp_tab[,8])
      }

      overview[[i]]$prior_prob     <- format(round(overview[[i]]$prior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$marg_lik       <- format(round(overview[[i]]$marg_lik, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$posterior_prob <- format(round(overview[[i]]$posterior_prob, x$add_info$digits_estimates), nsmall = x$add_info$digits_estimates)
      overview[[i]]$BF             <- format(round(overview[[i]]$BF, x$add_info$digits_BF), nsmall = x$add_info$digits_BF)
    }
  }


  cat("Call:\n")
  print(x$call)
  cat("\n")


  if(x$add_info$type == "ensemble"){

    cat("Robust Bayesian Survival Analysis\n")
    print(overview, quote = FALSE, right = TRUE)
    cat("\n")

    cat("Model-averaged estimates\n")
    print(averaged, quote = FALSE, right = TRUE)
    if(x$add_info$failed != 0)cat(paste0("\033[0;31m",x$add_info$failed, ifelse(x$add_info$failed == 1, " model", " models"), " failed to converge and ",ifelse(x$add_info$failed == 1, "was", "were")," omited from the summary.\033[0m\n"))


    if(!is.null(x$predictors$conditional)){
      cat("\n")
      cat("Conditional estimates\n")
      print(conditional, quote = FALSE, right = TRUE)
    }

    if(!is.null(x$intercept)){
      cat("\n")
      cat("Intercepts\n")
      print(intercept, quote = FALSE, right = TRUE)
    }

    if(!is.null(x$aux)){
      cat("\n")
      cat("Auxiliary parameters\n")
      print(aux, quote = FALSE, right = TRUE)
    }

  }else if(x$add_info$type == "models"){

    cat("Robust Bayesian Survival Analysis\n")
    print(overview, quote = FALSE, right = TRUE)

    if(!is.null(x$diagnostics)){
      cat("\n")
      cat("Models diagnostics overview\n")
      print(diagnostics, quote = FALSE, right = TRUE)
    }

  }else if(x$add_info$type == "individual"){

    cat("Individual Models Summary\n\n")

    for(i in 1:length(overview)){

      cat("Model Overview:\n")
      print(overview[[i]]$tab_main, quote = FALSE, right = TRUE, row.names = FALSE)

      cat("\nModel Coefficients:\n")
      print(overview[[i]]$tab_predictors, quote = FALSE, right = TRUE)

      cat("\nDistributional Parameters:\n")
      print(overview[[i]]$tab_int_and_aux, quote = FALSE, right = TRUE)
      if(i != length(overview))cat("\n\n")
    }
  }

}

#' @title Reports whether x is a RoBSA object
#'
#' @param x an object to test
#' @export is.RoBSA
is.RoBSA            <- function(x){
  inherits(x, "RoBSA")
}

