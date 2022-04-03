#' @title Check fitted RoBSA object for errors and warnings
#'
#' @description Checks fitted RoBSA object
#' for warnings and errors and prints them to the
#' console.
#'
#' @param fit a fitted RoBSA object.
#'
#'
#' @return \code{check_RoBSA} returns a vector of error and
#' warning messages.
#'
#' @export
check_RoBSA <- function(fit){
  .print_errors_and_warnings(fit, max_print = Inf)
}


.remove_model_posteriors   <- function(object){
  for(i in seq_along(object[["models"]])){
    if(inherits(object$models[[i]][["fit"]], "runjags")){
      object$models[[i]]$fit[["mcmc"]] <- NULL
    }
  }
  return(object)
}
.remove_model_margliks     <- function(object){
  for(i in seq_along(object[["models"]])){
    if(inherits(object$models[[i]][["marglik"]], "bridge")){
      object$models[[i]]$marglik[["q11"]] <- NULL
      object$models[[i]]$marglik[["q12"]] <- NULL
      object$models[[i]]$marglik[["q21"]] <- NULL
      object$models[[i]]$marglik[["q22"]] <- NULL
    }
  }
  return(object)
}
.print_errors_and_warnings <- function(object, max_print = 5){

  errors_and_warnings <- .collect_errors_and_warnings(object, max_print = max_print)

  for(i in seq_along(errors_and_warnings))
    warning(errors_and_warnings[i], immediate. = TRUE, call. = FALSE)

  return(invisible(errors_and_warnings))
}
.shorten_warnings    <- function(warnings, n_warnings = 5){
  if(is.null(warnings)){
    return(NULL)
  }else if(length(warnings) <= n_warnings){
    return(warnings)
  }else{
    return(c(warnings[1:n_warnings], paste0("There were another ", length(warnings) - n_warnings - 1, " warnings. To see all warnings call 'check_RoBSA(fit)'.")))
  }
}
.shorten_errors      <- function(errors, n_errors = 5){
  if(is.null(errors)){
    return(NULL)
  }else if(length(errors) <= n_errors){
    return(errors)
  }else{
    return(c(errors[1:n_errors], paste0("There were another ", length(errors) - n_errors - 1, " errors. To see all errors call 'check_RoBSA(fit)'.")))
  }
}
.convergence_warning <- function(object){
  if(any(!.get_model_convergence(object))){
    return(paste0(sum(!.get_model_convergence(object)), ifelse(sum(!.get_model_convergence(object)) == 1, " model", " models"), " failed to converge."))
  }else{
    return(NULL)
  }
}
.collect_errors_and_warnings <- function(object, max_print = 5){

  short_warnings <- .shorten_warnings(object$add_info[["warnings"]], max_print)
  short_errors   <- .shorten_errors(object$add_info[["errors"]],     max_print)
  conv_warning   <- .convergence_warning(object)

  return(c(short_warnings, short_errors, conv_warning))
}
.get_model_convergence       <- function(object){
  return(sapply(object[["models"]], function(model) if(is.null(model[["converged"]])) FALSE else model[["converged"]]))
}
.get_model_warnings          <- function(object){
  return(unlist(sapply(seq_along(object[["models"]]), function(i){
    if(!is.null(object[["models"]][[i]][["warnings"]])){
      paste0("Model (", i, "): ", object[["models"]][[i]][["warnings"]])
    }
  })))
}
.get_model_errors            <- function(object){
  return(unlist(sapply(seq_along(object[["models"]]), function(i){
    if(!is.null(object[["models"]][[i]][["errors"]])){
      paste0("Model (", i, "): ", object[["models"]][[i]][["errors"]])
    }
  })))
}


.has_aux         <- function(distributions){
  distributions != "exp-aft"
}
.intercept_name  <- function(distribution){
  switch(
    distribution,
    "exp-aft"     = "rate",
    "weibull-aft" = "scale",
    "lnorm-aft"   = "meanlog",
    "llogis-aft"  = "scale",
    "gamma-aft"   = "rate"
  )
}
.aux_name        <- function(distribution){
  switch(
    distribution,
    "exp-aft"     = NULL,
    "weibull-aft" = "shape",
    "lnorm-aft"   = "sdlog",
    "llogis-aft"  = "shape",
    "gamma-aft"   = "shape"
  )
}
.intercept_transformation <- function(distribution){
  switch(
    distribution,
    "exp-aft"     = function(x)exp(-x),
    "weibull-aft" = function(x)exp(x),
    "lnorm-aft"   = function(x)x,
    "llogis-aft"  = function(x)exp(x),
    "gamma-aft"   = function(x)exp(-x)
  )
}
.aux_transformation       <- function(distribution){
  switch(
    distribution,
    "exp-aft"     = NULL,
    "weibull-aft" = function(x)x,
    "lnorm-aft"   = function(x)x,
    "llogis-aft"  = function(x)x,
    "gamma-aft"   = function(x)x
  )
}
.is_zero_spike            <- function(prior){
  if(prior$distribution == "point"){
    if(prior$parameters$location == 0){
      return(TRUE)
    }
  }
  return(FALSE)
}
.is_model_constant        <- function(priors){
  # checks whether there is at least one non-nill prior
  return(all(sapply(priors, function(prior) BayesTools::is.prior.point(prior) | BayesTools::is.prior.none(prior))))
}

.reserved_words <- function() c("aux", "intercept", "terms", "time")
