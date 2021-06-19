# from RoBMA (with removal of some arguments)


#' @title Creates a RoBSA prior
#'
#' @description \code{prior} creates a prior distribution for fitting
#' a RoBSA model. The prior can be visualized by a \code{plot} function.
#'
#' @param distribution name of the prior distribution. The
#' possible options are
#' \describe{
#'   \item{\code{"point"}}{for a point density characterized by a
#'   \code{location} parameter.}
#'   \item{\code{"normal"}}{for a normal distribution characterized
#'   by a \code{mean} and \code{sd} parameters.}
#'   \item{\code{"lognormal"}}{for a lognormal distribution characterized
#'   by a \code{meanlog} and \code{sdlog} parameters.}
#'   \item{\code{"cauchy"}}{for a Cauchy distribution characterized
#'   by a \code{location} and \code{scale} parameters. Internally
#'   converted into a generalized t-distribution with \code{df = 1}.}
#'   \item{\code{"t"}}{for a generalized t-distribution characterized
#'   by a \code{location}, \code{scale}, and \code{df} parameters.}
#'   \item{\code{"gamma"}}{for a gamma distribution characterized
#'   by either \code{shape} and \code{rate}, or \code{shape} and
#'   \code{scale} parameters. The later is internally converted to
#'   the \code{shape} and \code{rate} parametrization}
#'   \item{\code{"invgamma"}}{for an inverse-gamma distribution
#'   characterized by a \code{shape} and \code{scale} parameters. The
#'   JAGS part uses a 1/gamma distribution with a shape and rate
#'   parameter.}
#'   \item{\code{"two.sided"}}{for a two-sided weight function
#'   characterized by a vector \code{steps} and vector \code{alpha}
#'   parameters. The \code{alpha} parameter determines an alpha
#'   parameter of Dirichlet distribution which cumulative sum
#'   is used for the weights omega.}
#'   \item{\code{"one.sided"}}{for a one-sided weight function
#'   characterized by either a vector \code{steps} and vector
#'   \code{alpha} parameter, leading to a monotonic one-sided
#'   function, or by a vector \code{steps}, vector \code{alpha1},
#'   and vector \code{alpha2} parameters leading non-monotonic
#'   one-sided weight function. The \code{alpha} / \code{alpha1} and
#'   \code{alpha2} parameters determine an alpha parameter of
#'   Dirichlet distribution which cumulative sum is used for
#'   the weights omega.}
#'   \item{\code{"uniform"}}{for a uniform distribution defined on a
#'   range from \code{a} to \code{b}}
#' }
#' @param parameters list of appropriate parameters for a given
#' \code{distribution}.
#' @param truncation list with two elements, \code{lower} and
#' \code{upper}, that define the lower and upper truncation of the
#' distribution. Defaults to \code{list(lower = -Inf, upper = Inf)}.
#' The lower truncation point is automatically set to 0 if it is
#' specified outside of the support of distributions defined only for
#' positive numbers.
#' @param prior_odds prior odds associated with a given distribution.
#' [RoBSA()] creates models corresponding to all combinations of prior
#' distributions for each of the model parameters (mu, tau, omega), and
#' sets the model priors odds to the product of its prior distributions.
#'
#' @examples
#' # create a standart normal prior distribution
#' p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))
#'
#' # create a half-normal standart normal prior distribution
#' p2 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1),
#' truncation = list(lower = 0, upper = Inf))
#'
#' # or a prior for one-sided weight function
#' p3 <- prior("one-sided", parameters = list(steps = c(.05, .10), alpha = c(1, 1, 1)))
#'
#' # the prior distribution can be visualized using the plot function
#' # (see ?plot.prior.RoBSA for all options)
#' plot(p1)
#'
#'
#' @export  prior
#' @rawNamespace S3method(print, RoBSA.prior)
#' @seealso [plot.RoBSA.prior()], \link[stats]{Normal}, \link[stats]{Cauchy},
#' \link[extraDistr]{LocationScaleT}, \link[stats]{GammaDist}, \link[extraDistr]{InvGamma}.
prior <- function(distribution, parameters, truncation = list(lower = -Inf, upper = Inf), prior_odds = 1){

  # general input check
  if(!(is.vector(distribution) & length(distribution) <= 2 & is.character(distribution)))stop("Argument 'distribution' is incorectly specified. It must be a character vector.")
  if(!is.list(parameters))stop("Argument 'parameters' must be a list.")
  if(!all(sapply(parameters, function(par)(is.vector(par) & is.numeric(par)))))stop("Elements of the 'parameter' argument must be numeric vectors.")
  if(!(is.list(truncation) & length(truncation) == 2))stop("Argument 'truncation' must be a list of length two.")
  if(!(is.vector(prior_odds) & is.numeric(prior_odds) & length(prior_odds) == 1))stop("Argument 'prior_odds' must be a numeric vector of length two.")
  if(is.null(names(truncation)))names(truncation) <- c("lower", "upper")
  if(truncation$lower >= truncation$upper)stop("The lower truncation point needs to be lower than the upper truncation point.")
  if(prior_odds <= 0)stop("Argument 'prior_odds' must be positive.")

  distribution <- tolower(distribution)

  if(length(distribution) == 2){
    if(any(distribution == "pet")){
      prefix       <- "PET"
      distribution <- distribution[distribution != "pet"]
    }else if(any(distribution == "peese")){
      prefix       <- "PEESE"
      distribution <- distribution[distribution != "peese"]
    }else{
      stop("The combination of distribution names was not recognized.")
    }
  }else{
    prefix <- NULL
  }

  # create an output object
  output <- list()

  # check whether the values are appropriate for the individual distribution, whether all parameters are included etc...
  if(distribution %in% c("norm", "normal")){

    if(length(parameters) != 2)stop("Normal prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("mean", "sd")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("mean", "sd")], sep = ", ", collapse = ""), " are not supported for a normal distribution."))
    }else{
      names(parameters) <- c("mean", "sd")
    }
    if(parameters$sd <= 0)stop("Parameter 'sd' must be positive.")

    # add the values to the output
    output$distribution <- "normal"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("lnorm", "lognormal")){

    if(length(parameters) != 2)stop("Lognormal prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("meanlog", "sdlog")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("meanlog", "sdlog")], sep = ", ", collapse = ""), " are not supported for a lognormal distribution."))
    }else{
      names(parameters) <- c("meanlog", "sdlog")
    }
    if(parameters$sdlog <= 0)stop("Parameter 'sdlog' must be positive.")

    if(truncation$lower == -Inf)truncation$lower <- 0 # change the defaul lower truncation
    if(truncation$lower < 0)stop("Lower truncation point must be non-negative.")

    # add the values to the output
    output$distribution <- "lognormal"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("t", "student", "student-t", "student t")){

    if(length(parameters) != 3)stop("Student-t prior distribution requires 3 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("df", "location", "scale")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("df", "location", "scale")], sep = ", ", collapse = ""), " are not supported for a student-t distribution."))
    }else{
      if(is.null(names(parameters)))names(parameters) <- c("location", "scale", "df")
    }
    if(parameters$df < 1)stop("Parameter 'df' must be positive.")
    if(!parameters$scale > 0)stop("Parameter 'scale' must be positive.")

    # add the values to the output
    output$distribution <- "t"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("cauchy")){

    if(length(parameters) != 2)stop("Cauchy prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("location", "scale")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("location", "scale")], sep = ", ", collapse = ""), " are not supported for a Cauchy distribution."))
    }else{
      if(is.null(names(parameters)))names(parameters) <- c("location", "scale")
    }
    if(!parameters$scale > 0)stop("Parameter 'scale' must be positive.")

    # add the values to the output
    output$distribution <- "t"
    output$parameters   <- list(location = parameters$location, scale = parameters$scale, df = 1) # pass as t-distribution
    output$truncation   <- truncation

  }else if(distribution %in% c("invgamma", "inversegamma", "inverse-gamma", "inverse gamma")){

    if(length(parameters) != 2)stop("Inverse gamma prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("shape", "scale")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("shape", "scale")], sep = ", ", collapse = ""), " are not supported for an inverse gamma distribution."))
    }else{
      names(parameters) <- c("shape", "scale")
    }
    if(!parameters$shape > 0)stop("Parameter 'shape' must be positive.")
    if(!parameters$scale > 0)stop("Parameter 'scale' must be positive.")

    if(truncation$lower == -Inf)truncation$lower <- 0 # change the defaul lower truncation
    if(truncation$lower < 0)stop("Lower truncation point must be non-negative.")
    # add the values to the output
    output$distribution <- "invgamma"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("gamma")){

    if(length(parameters) != 2)stop("Gamma prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("shape", "rate", "scale")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("shape", "rate", "scale")], sep = ", ", collapse = ""), " are not supported for a gamma distribution."))
    }else{
      names(parameters) <- c("shape", "rate")
    }
    if(!is.null(parameters$scale)){
      parameters$rate  <- 1/parameters$scale
      parameters$scale <- NULL
    }
    if(!parameters$shape > 0)stop("Parameter 'shape' must be positive.")
    if(!parameters$rate > 0)stop("Parameter 'rate' must be positive.")

    if(truncation$lower == -Inf)truncation$lower <- 0 # change the defaul lower truncation
    if(truncation$lower < 0)stop("Lower truncation point must be non-negative.")
    # add the values to the output
    output$distribution <- "gamma"
    output$parameters   <- parameters
    output$truncation   <- truncation

  }else if(distribution %in% c("point", "spike")){

    if(length(parameters) != 1)stop("Point prior distribution requires 1 parameter.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("location")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("location")], sep = ", ", collapse = ""), " are not supported for an inverse gamma distribution."))
    }else{
      names(parameters) <- c("location")
    }

    # add the values to the output
    output$distribution <- "point"
    output$parameters   <- parameters

  }else if(distribution %in% c("two.sided", "twosided", "two-sided", "two sided")){

    if(length(parameters) != 2)stop("Two-sided step function requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("alpha", "steps")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("alpha", "steps")], sep = ", ", collapse = ""), " are not supported for a two-sided weight function."))
    }else{
      names(parameters) <- c("steps", "alpha")
    }
    if(!all(parameters$alpha > 0))stop("Parameters 'alpha' must be positive.")
    if(!all(parameters$steps < 1 & parameters$steps > 0))stop("Parameters 'steps' must be higer than 0 and lower than 1.")
    if(!(all(parameters$steps == cummax(parameters$steps))))stop("Parameters 'steps' must be monotonically increasing.")
    if(length(parameters$steps) != length(parameters$alpha) - 1)stop("The parameter alpha needs to have one more argument then there are steps.")

    # reverse the ordering of alpha and weights - to corespond to ordering on t-statistics
    parameters$steps <- rev(parameters$steps)
    parameters$alpha <- rev(parameters$alpha)

    # add the values to the output
    output$distribution <- "two.sided"
    output$parameters   <- parameters

  }else if(distribution %in% c("one.sided", "onesided", "one-sided", "one sided")){

    if(!length(parameters) %in% c(2,3))stop("One-sided step function requires 2 or 3 parameters.")
    if(length(parameters) == 2){

      if(!is.null(names(parameters))){
        if(!all(names(parameters) %in% c("alpha", "steps")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("alpha", "steps")], sep = ", ", collapse = ""), " are not supported for a one-sided monotonic weight function."))
      }else{
        names(parameters) <- c("steps", "alpha")
      }
      if(!all(parameters$alpha > 0))stop("Parameters 'alpha' must be positive.")
      if(!all(parameters$steps < 1 & parameters$steps > 0))stop("Parameters 'steps' must be higer than 0 and lower than 1.")
      if(!(all(parameters$steps == cummax(parameters$steps))))stop("Parameters 'steps' must be monotonically increasing.")
      if(length(parameters$steps) != length(parameters$alpha) - 1)stop("The parameter alpha needs to have one more argument then there are steps.")

      # reverse the ordering of alpha and weights - to corespond to ordering on t-statistics
      parameters$steps <- rev(parameters$steps)
      parameters$alpha <- rev(parameters$alpha)

    }else if(length(parameters) == 3){

      if(!is.null(names(parameters))){
        if(!all(names(parameters) %in% c("alpha1", "alpha2", "steps")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("alpha1", "alpha2", "steps")], sep = ", ", collapse = ""), " are not supported for a one-sided non-monotonic weight function."))
      }else{
        names(parameters) <- c("steps","alpha1", "alpha2")
      }
      if(!all(parameters$alpha1 > 0))stop("Parameters 'alpha1' must be positive.")
      if(!all(parameters$alpha2 > 0))stop("Parameters 'alpha2' must be positive.")
      if(!all(parameters$steps < 1 & parameters$steps > 0))stop("Parameters 'steps' must be higer than 0 and lower than 1.")
      if(!(all(parameters$steps == cummax(parameters$steps))))stop("Parameters 'steps' must be monotonically increasing.")
      if(sum(parameters$steps <= .5) != length(parameters$alpha1) - 1)stop("The parameter alpha1 needs to have one more argument then there are steps <= .5.")
      if(sum(parameters$steps > .5) != length(parameters$alpha2) - 1)stop("The parameter alpha2 needs to have one more argument then there are steps > .5.")

      # reverse the ordering of alpha and weights - to corespond to ordering on t-statistics
      parameters$steps  <- rev(parameters$steps)
      parameters$alpha1 <- rev(parameters$alpha1)
    }

    # add the values to the output
    output$distribution <- "one.sided"
    output$parameters   <- parameters

  }else if(distribution %in% c("uniform", "unif")){

    if(length(parameters) != 2)stop("Uniform prior distribution requires 2 parameters.")
    if(!is.null(names(parameters))){
      if(!all(names(parameters) %in% c("a", "b")))stop(paste0("Parameters ", paste(names(parameters)[!names(parameters) %in% c("a", "b")], sep = ", ", collapse = ""), " are not supported for a normal distribution."))
    }else{
      names(parameters) <- c("a", "b")
    }
    if(parameters$a >= parameters$b)stop("Parameter 'a' must be lower than parameter 'b'.")
    if(!is.infinite(truncation$lower))if(truncation$lower != parameters$a)stop("Lower truncation must correspond to the parameter 'a'.")
    if(!is.infinite(truncation$upper))if(truncation$upper != parameters$b)stop("Upper truncation must correspond to the parameter 'b'.")

    # add the values to the output
    output$distribution <- "uniform"
    output$parameters   <- parameters
    output$truncation   <- list(lower = parameters$a, upper = parameters$b)

  }else{
    stop(paste0("The specified distribution name '", distribution,"' is not known. Please, see '?prior' for more information about supported prior distributions."))
  }

  # add a PET/PEESE prefix
  if(!is.null(prefix)){
    output$distribution <- paste0(prefix, ".", output$distribution)
  }

  output$prior_odds <- prior_odds
  class(output) <- "RoBSA.prior"

  return(output)
}


#' @title Prints a RoBSA.prior object
#'
#' @param x a RoBSA prior
#' @param ... additional arguments
#' \describe{
#'   \item{silent}{to silently return the print message.}
#'   \item{plot}{to return \link[base]{bquote} formatted
#'   prior name for plotting.}
#'   \item{digits_estimates}{number of decimals to be displayed
#'   for printed parameters.}
#'  }
#' @export  print.RoBSA.prior
#' @rawNamespace S3method(print, RoBSA.prior)
#' @seealso [prior()]
print.RoBSA.prior <- function(x, ...){

  dots <- list(...)
  silent            <- if(is.null(dots$silent))            FALSE else as.logical(dots$silent)
  plot              <- if(is.null(dots$plot))              FALSE else as.logical(dots$plot)
  digits_estimates  <- if(is.null(dots$digits_estimates))  2     else as.numeric(dots$digits_estimates)
  if(plot)silent <- TRUE


  # round the parameters and truncation for printing
  for(i in seq_along(x$parameters)){
    x$parameters[[i]] <- round(x$parameters[[i]], digits_estimates)
  }
  for(i in seq_along(x$truncation)){
    x$truncation[[i]] <- round(x$truncation[[i]], digits_estimates)
  }

  if(grepl("PET", x$distribution)){
    prefix          <- "PET:"
    x$distribution  <- substr(x$distribution, 5, nchar(x$distribution))
  }else if(grepl("PEESE", x$distribution)){
    prefix          <- "PEESE:"
    x$distribution  <- substr(x$distribution, 7, nchar(x$distribution))
  }else{
    prefix <- NULL
  }

  name <- switch(x$distribution,
                 "normal"       = "Normal",
                 "lognormal"    = "LogNormal",
                 "t"            = "gen. Student-t",
                 "gamma"        = "Gamma",
                 "invgamma"     = "InvGamma",
                 "point"        = "Spike",
                 "two.sided"    = "Two-sided",
                 "one.sided"    = "One-sided",
                 "uniform"      = "Uniform")

  name <- paste0(prefix, name)

  if(x$distribution %in% c("normal", "PET.normal", "PEESE.normal")){
    parameters <- c(x$parameters$mean, x$parameters$sd)
  }else if(x$distribution == "lognormal"){
    parameters <- c(x$parameters$meanlog, x$parameters$sdlog)
  }else if(x$distribution == "t"){
    if(x$parameters$df == 1){
      name <- "Cauchy"
      parameters <- c(x$parameters$location, x$parameters$scale)
    }else{
      parameters <- c(x$parameters$location, x$parameters$scale, x$parameters$df)
    }
  }else if(x$distribution == "gamma"){
    parameters <- c(x$parameters$shape, x$parameters$rate)
  }else if(x$distribution == "invgamma"){
    parameters <- c(x$parameters$shape, x$parameters$scale)
  }else if(x$distribution == "point"){
    parameters <- c(x$parameters$location)
  }else if(x$distribution == "two.sided"){
    parameters <- c(
      paste0("(",paste(x$parameters$steps, collapse = ", "),")"),
      paste0("(",paste(x$parameters$alpha, collapse = ", "),")")
    )
  }else if(x$distribution == "one.sided"){
    if(all(names(x$parameters) %in% c("steps", "alpha1", "alpha2"))){
      parameters <- c(
        paste0("(",paste(x$parameters$steps, collapse = ", "),")"),
        paste0("(",paste(x$parameters$alpha1, collapse = ", "),")"),
        paste0("(",paste(x$parameters$alpha2, collapse = ", "),")")
      )
    }else{
      parameters <- c(
        paste0("(",paste(x$parameters$steps, collapse = ", "),")"),
        paste0("(",paste(x$parameters$alpha, collapse = ", "),")")
      )
    }

  }else if(x$distribution == "uniform"){
    parameters <- c(x$parameters$a, x$parameters$b)
  }

  if(plot){
    if(!x$distribution %in% c("point","one.sided","two.sided","uniform")){

      output <- bquote(italic(.(name))*.(paste0("(",paste(parameters, collapse = ", "),")"))[~"["~
                                                                                               .(if(is.infinite(x$truncation$lower)){bquote(-infinity)}else{x$truncation$lower})*", "*
                                                                                               .(if(is.infinite(x$truncation$upper)){bquote( infinity)}else{x$truncation$upper})~"]"])
    }else{
      output <- bquote(italic(.(name))*.(paste0("(",paste(parameters, collapse = ", "),")")))
    }

  }else{
    output <- paste0(name,"(",paste(parameters, collapse = ", "),")")
    if(!x$distribution %in% c("point","one.sided","two.sided","uniform")){
      output <- paste0(output, "[",x$truncation$lower, ", ",x$truncation$upper, "]")
    }
  }


  if(!silent)cat(output)
  return(invisible(output))
}

