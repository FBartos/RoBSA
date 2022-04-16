#' @title Plots a fitted RoBSA object
#'
#' @description \code{plot.RoBSA} allows to visualize
#' posterior distribution of different \code{"RoBSA"} object
#' parameters. See \code{plot_survival} for plotting the survival
#' ways. See \code{type} for the different model types.
#'
#' @param x a fitted RoBSA object
#' @param parameter a name of parameter to be plotted. Defaults to
#' the first regression parameter if left unspecified. Use
#' \code{"intercept"} and \code{"aux"} to plot the intercepts and
#' auxiliary parameters of each distribution family.
#' @param conditional whether conditional estimates should be
#' plotted. Defaults to \code{FALSE} which plots the model-averaged
#' estimates.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # the 'plot' function allows to visualize the results of a fitted RoBSA object, for example;
#' # the model-averaged effect size estimate
#' }
#'
#'
#' @return \code{plot.RoBSA} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBSA()]
#' @export
plot.RoBSA  <- function(x, parameter = NULL, conditional = FALSE, plot_type = "base", prior = FALSE, dots_prior = NULL, ...){

  # check whether plotting is possible
  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")

  #check settings
  BayesTools::check_char(parameter, "parameter", allow_NULL = TRUE)
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")

  # verify whether given parameter exists and set a default parameter for plotting
  if(is.null(parameter)){
    parameter <- x$add_info[["predictors"]][1]
  }else if(!parameter %in% c("intercept", "aux", x$add_info[["predictors"]])){
    stop(paste0("The passed parameter does not correspond to any of the specified predictors: ", paste0("'", x$add_info[["predictors"]], "'", collapse = ", ")))
  }


  # choose the samples
  if(parameter == "intercept"){
    samples <- x[["RoBSA"]][["posteriors_intercept"]]
  }else if(parameter == "aux"){
    samples <- x[["RoBSA"]][["posteriors_aux"]]
  }else{
    if(conditional){
      samples <- x[["RoBSA"]][["posteriors_conditional"]]
    }else{
      samples <- x[["RoBSA"]][["posteriors"]]
    }
  }


  dots       <- .set_dots_plot(...)
  dots_prior <- .set_dots_prior(dots_prior)

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$transformation           <- NULL
  args$transformation_arguments <- NULL
  args$transformation_settings  <- FALSE
  args$rescale_x                <- FALSE
  args$dots_prior               <- dots_prior

  if(parameter %in% c("intercept", "aux")){
    plot <- list()
    for(i in seq_along(samples)){
      args$parameter <- names(samples)[i]
      args$par_name  <- names(samples)[i]
      plot[[names(samples)[i]]] <- do.call(BayesTools::plot_posterior, args)
    }
  }else{
    args$parameter <- .BayesTools_parameter_name(parameter)
    args$par_name  <- parameter
    plot           <- do.call(BayesTools::plot_posterior, args)
  }


  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


.set_dots_plot        <- function(...){

  dots <- list(...)
  if(is.null(dots[["col"]])){
    dots[["col"]]      <- "black"
  }
  if(is.null(dots[["col.fill"]])){
    dots[["col.fill"]] <- "#4D4D4D4C" # scales::alpha("grey30", .30)
  }

  return(dots)
}
.set_dots_prior       <- function(dots_prior){

  if(is.null(dots_prior)){
    dots_prior <- list()
  }

  if(is.null(dots_prior[["col"]])){
    dots_prior[["col"]]      <- "grey60"
  }
  if(is.null(dots_prior[["lty"]])){
    dots_prior[["lty"]]      <- 1
  }
  if(is.null(dots_prior[["col.fill"]])){
    dots_prior[["col.fill"]] <- "#B3B3B34C" # scales::alpha("grey70", .30)
  }

  return(dots_prior)
}


#' @title Models plot for a RoBSA object
#'
#' @description \code{plot_models} plots individual models'
#' estimates for a \code{"RoBSA"} object.
#'
#' @param parameter a name of parameter to be plotted. Defaults to
#' the first regression parameter if left unspecified.
#' @param order how the models should be ordered.
#' Defaults to \code{"decreasing"} which orders them in decreasing
#' order in accordance to \code{order_by} argument. The alternative is
#' \code{"increasing"}.
#' @param order_by what feature should be use to order the models.
#' Defaults to \code{"model"} which orders the models according to
#' their number. The alternatives are \code{"estimate"} (for the effect
#' size estimates), \code{"probability"} (for the posterior model probability),
#' and \code{"BF"} (for the inclusion Bayes factor).
#' @inheritParams plot.RoBSA
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBSA(r = Anderson2010$r, n = Anderson2010$n, study_names = Anderson2010$labels)
#'
#' ### ggplot2 version of all of the plots can be obtained by adding 'model_type = "ggplot"
#' # the plot_models function creates a plot for of the individual models' estimates, for example,
#' # the effect size estimates from the individual models can be obtained with
#' plot_models(fit)
#'
#' # and effect size estimates from only the conditional models
#' plot_models(fit, conditional = TRUE)
#' }
#'
#'
#' @return \code{plot_models} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @export
plot_models <- function(x, parameter = NULL, conditional = FALSE, plot_type = "base", order = "decreasing", order_by = "model", ...){

  if(sum(.get_model_convergence(x)) == 0)
    stop("There is no converged model in the ensemble.")

  #check settings
  BayesTools::check_char(parameter, "parameter", allow_NULL = TRUE)
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_char(order, "order", allow_NULL = TRUE, allow_values = c("increasing", "decreasing"))
  BayesTools::check_char(order, "order", allow_NULL = TRUE, allow_values = c("increasing", "decreasing"))
  BayesTools::check_char(order_by, "order_by", allow_NULL = TRUE, allow_values = c("model", "estimate", "probability", "BF"))


  # verify whether given parameter exists and set a default parameter for plotting
  if(is.null(parameter)){
    parameter <- x$add_info[["predictors"]][1]
  }else if(!parameter %in% x$add_info[["predictors"]]){
    stop(paste0("The passed parameter does not correspond to any of the specified predictors: ", paste0("'", x$add_info[["predictors"]], "'", collapse = ", ")))
  }


  ### prepare input
  if(conditional){

    model_list <- x[["models"]]
    samples    <- x[["RoBSA"]][["posteriors_conditional"]]
    inference  <- x[["RoBSA"]][["inference_conditional"]]

  }else{

    model_list <- x[["models"]]
    samples    <- x[["RoBSA"]][["posteriors"]]
    inference  <- x[["RoBSA"]][["inference"]]

  }

  dots <- list(...)

  # prepare the argument call
  args                          <- dots
  args$model_list               <- model_list
  args$samples                  <- samples
  args$inference                <- inference
  args$parameter                <- .BayesTools_parameter_name(parameter)
  args$par_name                 <- NULL
  args$plot_type                <- plot_type
  args$prior                    <- FALSE
  args$conditional              <- conditional
  args$order                    <- list(order, order_by)
  args$transformation           <- NULL
  args$transformation_arguments <- NULL
  args$transformation_settings  <- FALSE
  args$formula_prefix           <- FALSE

  plot <- do.call(BayesTools::plot_models, args)


  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


#' @title Plots a fitted RoBSA object
#'
#' @param x a fitted RoBSA object.
#' @param ... additional arguments.
#' @export
plot_survival <- function(x, parameter, conditional, data_covariates = NULL, time_range = NULL){

  # create times for the prediction
  df          <- x$data
  time_range  <- c(0, max(c(df$t_event, df$tcens_r)))
  pred_data   <- data.frame(time = seq(time_range[1], time_range[2], length.out = 100))

  # add mean values for the predictors
  pred_names  <- attr(df, "predictors")
  for(i in seq_along(pred_names)){
    mean_pred   <- mean(unlist(df[grep(pred_names[i], names(df))]))
    pred_data   <- cbind(pred_data, new_pred = rep(mean_pred, nrow(pred_data)))
    colnames(pred_data)[ncol(pred_data)] <- pred_names[i]
  }

  # predict the mean survival
  pred_survival <- predict.RoBSA(x, pred_data, type = "survival")

  # plot the survival
  graphics::plot(x = pred_data$time, y = pred_survival$mean, ylim = c(0,1), "l", xlab = "Time", ylab = "Survival", las = 1, ...)
  graphics::lines(x = pred_data$time, y = pred_survival$uCI, lty = 2, ...)
  graphics::lines(x = pred_data$time, y = pred_survival$lCI, lty = 2, ...)

  return(invisible())
}
