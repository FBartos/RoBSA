#' @title Visualizes MCMC diagnostics for a fitted RoBSA object
#'
#' @description \code{diagnostics} creates visual
#' checks of individual models convergence. Numerical
#' overview of individual models can be obtained by
#' \code{summary(object, type = "diagnostics")},
#' or even more detailed information by
#' \code{summary(object, type = "individual")}.
#'
#' @param fit a fitted RoBSA object
#' @param parameter a parameter to be plotted.
#' @param type type of MCMC diagnostic to be plotted.
#' Options are \code{"trace"} for the chains' trace plots,
#' \code{"autocorrelation"} for autocorrelation of the
#' chains, and \code{"densities"} for the overlaying
#' densities of the individual chains. Can be abbreviated to
#' first letters.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param show_models MCMC diagnostics of which models should be
#' plotted. Defaults to \code{NULL} which plots MCMC diagnostics
#' for a specified parameter for every model that is part of the
#' ensemble.
#' @param title whether the model number should be displayed in title.
#' Defaults to \code{TRUE} when more than one model is selected.
#' @param lags number of lags to be shown for
#' \code{type = "autocorrelation"}. Defaults to \code{30}.
#' @param ... additional arguments to be passed to the plotting functions.
#'
#'
#' @return \code{diagnostics} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object/list of objects (depending on the number of parameters to be plotted)
#' of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBSA()], [summary.RoBSA()]
#'
#' @name diagnostics
#' @aliases diagnostics_autocorrelation diagnostics_trace diagnostics_density
#' @export diagnostics
#' @export diagnostics_density
#' @export diagnostics_autocorrelation
#' @export diagnostics_trace

#' @rdname diagnostics
diagnostics <- function(fit, parameter = NULL, type, plot_type = "base", show_models = NULL,
                        lags = 30, title = is.null(show_models) | length(show_models) > 1, ...){

  # check settings
  if(class(fit) != "RoBSA")
    stop("Diagnostics are available only for RoBSA models.")
  if(fit$add_info[["save"]] == "min")
    stop("Diagnostics cannot be produced because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' while fitting the model (see ?RoBSA for more details).")
  BayesTools::check_char(parameter, "parameter", allow_NULL = TRUE)
  BayesTools::check_char(type, "type", allow_values = c("density", "trace", "autocorrelation"))
  BayesTools::check_bool(title, "title", allow_NULL = TRUE)

  # verify whether given parameter exists and set a default parameter for plotting
  if(is.null(parameter)){
    parameter <- fit$add_info[["predictors"]][1]
  }else if(!parameter %in% c("intercept", "aux", fit$add_info[["predictors"]])){
    stop(paste0("The passed parameter does not correspond to any of the specified predictors: ", paste0("'", fit$add_info[["predictors"]], "'", collapse = ", ")))
  }


  # obtain model indexes to plot
  models_ind <- 1:length(fit[["models"]])
  if(!is.null(show_models)){
    models_ind <- models_ind[show_models]
  }

  # a message with info about multiple plots
  if(plot_type == "base" & length(models_ind) > 1)
    message("Multiple plots will be produced. See '?layout' for help with setting multiple plots.")


  dots  <- .set_dots_diagnostics(..., type = type, chains = fit[["fit_control"]][["chains"]])
  plots <- list()

  for(i in models_ind){

    if(!parameter %in% c("intercept", "aux", fit$models[[i]][["terms_test"]])){

      plots[[i]] <- NULL

    }else{

      # get the parameter name
      args                          <- dots
      args$fit                      <- fit$models[[i]][["fit"]]
      args$parameter                <- if(parameter == "aux") "aux" else .BayesTools_parameter_name(parameter)
      args$type                     <- type
      args$plot_type                <- plot_type
      args$lags                     <- lags
      args$transformations          <- NULL
      args$transform_orthonormal    <- TRUE
      args$short_name               <- FALSE
      args$parameter_names          <- FALSE
      args$formula_prefix           <- FALSE

      if(!is.null(title) && title){
        args$main <- paste0("Model ", i)
      }

      plots[[i]] <- do.call(BayesTools::JAGS_diagnostics, args)
    }
  }



  # return the plots
  if(plot_type == "base"){
    return(invisible(plots))
  }else if(plot_type == "ggplot"){
    if(length(plots) == 1){
      plots <- plots[[1]]
    }
    return(plots)
  }
}

#' @rdname diagnostics
diagnostics_autocorrelation <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        lags = 30, title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "autocorrelation", plot_type = plot_type, show_models = show_models, lags = lags, title = title, ...)
}

#' @rdname diagnostics
diagnostics_trace           <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "trace", plot_type = plot_type, show_models = show_models, title = title, ...)
}

#' @rdname diagnostics
diagnostics_density         <- function(fit, parameter = NULL, plot_type = "base", show_models = NULL,
                                        title = is.null(show_models) | length(show_models) > 1, ...){
  diagnostics(fit = fit, parameter = parameter, type = "density", plot_type = plot_type, show_models = show_models, title = title, ...)
}

.set_dots_diagnostics  <- function(..., type, chains){

  dots <- list(...)
  if(is.null(dots[["col"]])){
    dots[["col"]]      <- if(type == "autocorrelation") "black" else rev(scales::viridis_pal()(chains))
  }

  return(dots)
}
