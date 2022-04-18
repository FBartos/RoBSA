#' @title Extract parameter estimates from \code{flexsurv} object
#'
#' @description \code{extract_flexsurv} extracts estimates from a
#' \code{flexsurv} object in and transform them so they match the
#' \code{RoBSA} output.
#'
#' @param fit an object fitted with the \code{flexsurv::flexsurvreg}
#' function
#'
#' @return \code{extract_flexsurv} return list of estimates lists for each
#' parameter.
#'
#' @export
extract_flexsurv <- function(fit){

  if(!inherits(fit, "flexsurvreg"))
    stop("'fit' must be a model fitted with the 'flexsurv::flexsurvreg' function.")

  distribution <- .flexsurv_distribution_name(fit)
  estimates    <- fit$res

  out <- list()


  out[["intercept"]] <- list(
      "mean" = .flexsurv_intercept_transform(distribution)(estimates[.flexsurv_intercept_name(distribution), "est"]),
      "lCI"  = .flexsurv_intercept_transform(distribution)(estimates[.flexsurv_intercept_name(distribution), "est"] -
                                                             1.96 * estimates[.flexsurv_intercept_name(distribution), "se"]),
      "uCI"  = .flexsurv_intercept_transform(distribution)(estimates[.flexsurv_intercept_name(distribution), "est"] +
                                                             1.96 * estimates[.flexsurv_intercept_name(distribution), "se"])
  )

  if(.has_aux(distribution)){
    out[["aux"]] <- list(
      "mean" = .flexsurv_aux_transform(distribution)(estimates[.flexsurv_aux_name(distribution), "est"]),
      "lCI"  = .flexsurv_aux_transform(distribution)(estimates[.flexsurv_aux_name(distribution), "est"] -
                                                       1.96 * estimates[.flexsurv_aux_name(distribution), "se"]),
      "uCI"  = .flexsurv_aux_transform(distribution)(estimates[.flexsurv_aux_name(distribution), "est"] +
                                                       1.96 * estimates[.flexsurv_aux_name(distribution), "se"])
    )
  }

  # omit parameter names
  estimates <- estimates[!rownames(estimates) %in% c(.flexsurv_intercept_name(distribution), .flexsurv_aux_name(distribution)),,drop=FALSE]
  if(nrow(estimates) > 0){
    out[["parameters"]]        <- lapply(rownames(estimates), function(est) list(
      "mean" = .flexsurv_est_transform(distribution)(estimates[est, "est"]),
      "lCI"  = .flexsurv_est_transform(distribution)(estimates[est, "est"] - 1.96 * estimates[est, "se"]),
      "uCI"  = .flexsurv_est_transform(distribution)(estimates[est, "est"] + 1.96 * estimates[est, "se"])
    ))
    names(out[["parameters"]]) <- rownames(estimates)
  }

  return(out)
}

.flexsurv_distribution_name    <- function(fit){
  switch(
    fit$dlist$name,
    "weibull.quiet" = "weibull-aft",
    "lnorm"         = "lnorm-aft",
    "llogis"        = "llogis-aft",
    "exp"           = "exp-aft",
    "gamma"         = "gamma-aft"
  )
}
.flexsurv_intercept_name       <- function(distribution){
  switch(
    distribution,
    "weibull-aft"   = "scale",
    "lnorm-aft"     = "meanlog",
    "llogis-aft"    = "scale",
    "exp-aft"       = "rate",
    "gamma-aft"     = "rate"
  )
}
.flexsurv_intercept_transform  <- function(distribution){
  switch(
    distribution,
    "weibull-aft"   = function(x) log(x),
    "lnorm-aft"     = function(x) x,
    "llogis-aft"    = function(x) log(x),
    "exp-aft"       = function(x) log(1/x),
    "gamma-aft"     = function(x) x
  )
}
.flexsurv_aux_name             <- function(distribution){
  switch(
    distribution,
    "weibull-aft"   = "shape",
    "lnorm-aft"     = "sdlog",
    "llogis-aft"    = "shape",
    "gamma-aft"     = "shape"
  )
}
.flexsurv_aux_transform        <- function(distribution){
  switch(
    distribution,
    "weibull-aft"   = function(x) x,
    "lnorm-aft"     = function(x) x,
    "llogis-aft"    = function(x) x,
    "gamma-aft"     = function(x) x
  )
}
.flexsurv_est_transform        <- function(distribution){
  switch(
    distribution,
    "exp-aft"       = function(x) -x,
    "weibull-aft"   = function(x) x,
    "lnorm-aft"     = function(x) x,
    "llogis-aft"    = function(x) x,
    "gamma-aft"     = function(x) -x
  )
}
