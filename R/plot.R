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
#' # plot posterior distribution of the treatment effect
#' plot(fit, parameter = "treatment")
#'
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

  # apply version changes to RoBSA object
  x <- .update_object(x)

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
#' # plot posterior distribution of the treatment effect from each model
#' plot_models(fit, parameter = "treatment")
#'
#' }
#'
#'
#' @return \code{plot_models} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @export
plot_models <- function(x, parameter = NULL, conditional = FALSE, plot_type = "base", order = "decreasing", order_by = "model", ...){

  if(!is.RoBSA(x))
    stop("'plot_models' function supports only 'RoBSA' type objects.")
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


#' @title Survival plots for a RoBSA object
#'
#' @param x a fitted RoBSA object.
#' @param time_range a numeric of length two specifying the range for the
#' survival prediction. Defaults to \code{NULL} which uses the range of
#' observed times.
#'
#' @inheritParams predict.RoBSA
#' @inheritParams plot.RoBSA
#' @param ... additional arguments.
#'
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
#' # plot survival for each level the treatment
#' plot_survival(fit, parameter = "treatment")
#'
#' # plot hazard for each level the treatment
#' plot_hazard(fit, parameter = "treatment")
#'
#' # plot density for each level the treatment
#' plot_density(fit, parameter = "treatment")
#' }
#'
#' @return returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @export plot_prediction
#' @export plot_survival
#' @export plot_hazard
#' @export plot_density
#' @name plot_prediction
#' @aliases plot_survival plot_hazard plot_density
NULL

#' @rdname plot_prediction
plot_prediction <- function(x, type = "survival", time_range = NULL, new_data = NULL, predictor = NULL, covariates_data = NULL,
                            conditional = FALSE, plot_type = "base", samples = 10000, ...){

  if(!is.RoBSA(x))
    stop("'plot_models' function supports only 'RoBSA' type objects.")
  BayesTools::check_real(time_range, "time_range", lower = 0, allow_NULL = TRUE, check_length = TRUE)
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  # the other input is checked within predict.RoBSA

  # set a default time_range
  if(is.null(time_range)){
    time_range <- c(0, max(c(
      x[["data"]][["survival"]][["t_event"]],
      x[["data"]][["survival"]][["t_cens_r"]],
      x[["data"]][["survival"]][["t_cens_ir"]]
    )))
    time_range <- range(pretty(time_range))
  }
  time_range <- seq(time_range[1], time_range[2], length.out = 101)
  if(type == "density"){
    time_range <- time_range[-1]
  }

  # obtain predictions for the survival
  plot_data  <- predict.RoBSA(x, time = time_range, new_data = new_data, predictor = predictor, covariates_data = covariates_data,
                              type = type, summarize = TRUE, averaged = TRUE, conditional = conditional, samples = samples)

  # generate plotting settings
  dots <- .set_dots_plot_prediction(plot_data, ...)


  # create times for the prediction
  if(plot_type == "base"){

    plot <- NULL
    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = dots[["xlab"]], ylab = dots[["ylab"]], main = dots[["main"]],
                   xlim = dots[["xlim"]], ylim = dots[["ylim"]],
                   cex.axis = dots[["cex.axis"]], cex.lab = dots[["cex.lab"]], cex.main = dots[["cex.main"]],
                   col.axis = dots[["col.axis"]], col.lab = dots[["col.lab"]], col.main = dots[["col.main"]])

    for(i in seq_along(plot_data)){
      graphics::polygon(
        x   = c(plot_data[[i]][,"time"], rev(plot_data[[i]][,"time"])),
        y   = c(plot_data[[i]][,"lCI"],  rev(plot_data[[i]][,"uCI"])),
        col = dots[["col.fill"]][i], border = NA
      )
    }

    for(i in seq_along(plot_data)){
      graphics::lines(plot_data[[i]][,"time"], plot_data[[i]][,"mean"], lwd = dots[["lwd"]][i], lty = dots[["lty"]][i], col = dots[["col"]][i])
    }

    if(dots[["legend"]]){
      graphics::legend(
        dots[["legend.position"]],
        legend = dots[["legend.text"]],
        col    = dots[["col"]],
        lty    = dots[["lty"]],
        lwd    = dots[["lwd"]],
        bty    = "n")
    }

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()+
      ggplot2::ggtitle(dots[["main"]]) +
      ggplot2::scale_x_continuous(name = dots[["xlab"]], limits = dots[["xlim"]], oob = scales::oob_keep) +
      ggplot2::scale_y_continuous(name = dots[["ylab"]], limits = dots[["ylim"]], oob = scales::oob_keep)

    for(i in seq_along(plot_data)){
      plot <- plot + ggplot2::geom_polygon(
        data    = data.frame(
          x = c(plot_data[[i]][,"time"], rev(plot_data[[i]][,"time"])),
          y = c(plot_data[[i]][,"lCI"],  rev(plot_data[[i]][,"uCI"]))),
        mapping = ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]]),
        fill    = dots[["col.fill"]][i]
      )
    }

    for(i in seq_along(plot_data)){
      plot <- plot + ggplot2::geom_line(
        data    = data.frame(
          x     = plot_data[[i]][,"time"],
          y     = plot_data[[i]][,"mean"],
          level = if(dots[["legend"]]) dots[["legend.text"]][i] else ""),
        mapping = ggplot2::aes(
          x         = .data[["x"]],
          y         = .data[["y"]],
          color     = if(length(dots[["legend.text"]]) > 1) .data[["level"]],
          linetype  = if(length(dots[["legend.text"]]) > 1) .data[["level"]],
          group     = if(length(dots[["legend.text"]]) > 1) .data[["level"]]),
        show.legend = dots[["legend"]])
    }

    if(dots[["legend"]]){
      names(dots[["lty"]]) <- dots[["legend.text"]]
      names(dots[["col"]]) <- dots[["legend.text"]]
      names(dots[["lwd"]]) <- dots[["legend.text"]]
      plot <- plot +
        ggplot2::scale_linetype_manual(name = "level", values = dots[["lty"]]) +
        ggplot2::scale_color_manual(name = "level", values = dots[["col"]]) +
        ggplot2::scale_size_manual(name = "level", values = dots[["lwd"]]) +
        ggplot2::theme(
        legend.title    = ggplot2::element_blank(),
        legend.position = dots[["legend_position"]])
    }
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
#' @rdname plot_prediction
plot_survival   <- function(x, time_range = NULL, new_data = NULL, predictor = NULL, covariates_data = NULL,
                            conditional = FALSE, plot_type = "base", samples = 10000, ...){
  plot_prediction(x, type = "survival", time_range = time_range, new_data = new_data, predictor = predictor, covariates_data = covariates_data,
                  conditional = conditional, plot_type = plot_type, samples = samples, ...)
}
#' @rdname plot_prediction
plot_hazard     <- function(x, time_range = NULL, new_data = NULL, predictor = NULL, covariates_data = NULL,
                            conditional = FALSE, plot_type = "base", samples = 10000, ...){
  plot_prediction(x, type = "hazard", time_range = time_range, new_data = new_data, predictor = predictor, covariates_data = covariates_data,
                  conditional = conditional, plot_type = plot_type, samples = samples, ...)
}
#' @rdname plot_prediction
plot_density    <- function(x, time_range = NULL, new_data = NULL, predictor = NULL, covariates_data = NULL,
                            conditional = FALSE, plot_type = "base", samples = 10000, ...){
  plot_prediction(x, type = "density", time_range = time_range, new_data = new_data, predictor = predictor, covariates_data = covariates_data,
                  conditional = conditional, plot_type = plot_type, samples = samples, ...)
}

.set_dots_plot_prediction <- function(plot_data, outcome, ...){

  dots   <- list(...)
  data   <- attr(plot_data, "data")
  levels <- nrow(data)

  if(is.null(dots[["xlab"]])){
    dots[["xlab"]] <- "Time"
  }

  if(is.null(dots[["ylab"]])){
    dots[["ylab"]] <- switch(
      attr(plot_data, "outcome"),
      "survival" = "Survival",
      "hazard"   = "Hazard",
      "density"  = "Density"
    )
  }

  if(is.null(dots[["main"]])){
    dots[["main"]] <- NULL
  }

  if(is.null(dots[["xlim"]])){
    dots[["xlim"]] <- range(attr(plot_data, "time"))
  }

  if(is.null(dots[["ylim"]])){
    dots[["ylim"]] <- switch(
      attr(plot_data, "outcome"),
      "survival" = c(0, 1),
      "hazard"   = c(0, 1),
      "density"  = range(pretty(range(do.call(c, lapply(plot_data, function(d) d[,-1])))))
    )
  }

  if(is.null(dots[["col"]])){
    if(levels == 1){
      dots[["col"]] <- "black"
    }else{
      dots[["col"]] <- scales::viridis_pal()(levels)
    }
  }else if(length(dots[["col"]]) == 1){
    dots[["col"]] <- rep(dots[["col"]], levels)
  }else if(dots[["col"]] != levels){
    stop("The number of specified colors does not match the number of predicted variable levels.")
  }

  if(is.null(dots[["col.fill"]])){
    if(levels == 1){
      dots[["col.fill"]] <- scales::alpha("grey30", alpha = .30)
    }else{
      dots[["col.fill"]] <- scales::alpha(dots[["col"]], alpha = .30)
    }
  }else if(length(dots[["col.fill"]]) == 1){
    dots[["col.fill"]] <- rep(dots[["col.fill"]], levels)
  }else if(dots[["col.fill"]] != levels){
    stop("The number of specified filling colors does not match the number of predicted variable levels.")
  }

  if(is.null(dots[["lwd"]])){
    dots[["lwd"]] <- rep(1, levels)
  }else if(length(dots[["lwd"]]) == 1){
    dots[["lwd"]] <- rep(dots[["lwd"]], levels)
  }else if(dots[["lwd"]] != levels){
    stop("The number of specified line widths (lwd) colors does not match the number of predicted variable levels.")
  }

  if(is.null(dots[["lty"]])){
    dots[["lty"]] <- rep(1, levels)
  }else if(length(dots[["lty"]]) == 1){
    dots[["lty"]] <- rep(dots[["lty"]], levels)
  }else if(dots[["lty"]] != levels){
    stop("The number of specified line types (lty) colors does not match the number of predicted variable levels.")
  }

  if(is.null(dots[["legend"]])){
    dots[["legend"]] <- TRUE
  }

  if(is.null(dots[["legend.text"]])){
    # remove constant predictors for defaults legend
    for(i in ncol(data):1){
      if(length(unique(data[,i])) == 1){
        data <- data[,-i, drop = FALSE]
      }
    }
    if(ncol(data) >= 1){
      dots[["legend.text"]] <- sapply(1:nrow(data), function(i) paste0(
        colnames(data), " = ", sapply(data[i,], function(x) if(is.numeric(x)) sprintf("%.2f", x) else x),
        collapse = "; "))
    }
  }

  if(is.null(dots[["legend.position"]])){
    dots[["legend.position"]] <- "topright"
  }

  if(length(dots[["legend.text"]]) == 0){
    dots[["legend"]] <- FALSE
  }


  return(dots)
}
