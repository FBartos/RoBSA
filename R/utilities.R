#' @title Options for the RoBSA package
#'
#' @description A placeholder object and functions for the RoBSA package.
#' (adapted from the runjags R package).
#'
#' @param name the name of the option to get the current value of - for a list of
#' available options, see details below.
#' @param ... named option(s) to change - for a list of available options, see
#' details below.
#'
#' @return The current value of all available RoBSA options (after applying any
#' changes specified) is returned invisibly as a named list.
#'
#' @export RoBSA.options
#' @export RoBSA.get_option
#' @name RoBSA_options
#' @aliases RoBSA_options RoBSA.options RoBSA.get_option
NULL


#' @rdname RoBSA_options
RoBSA.options    <- function(...){

	opts <- list(...)

	for(i in seq_along(opts)){

	  if(!names(opts)[i] %in% names(RoBSA.private))
	    stop(paste("Unmatched or ambiguous option '", names(opts)[i], "'", sep=""))

	  if(names(opts)[i] == "distribution")
	    stop("List of default distributions cannot be changed.")

	  assign(names(opts)[i], opts[[i]] , envir = RoBSA.private)
	}

	return(invisible(RoBSA.private$options))
}

#' @rdname RoBSA_options
RoBSA.get_option <- function(name){

	if(length(name)!=1)
	  stop("Only 1 option can be retrieved at a time")

	if(!name %in% names(RoBSA.private))
	  stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))

	# Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
	return(eval(RoBSA.private[[name]]))
}

# export the function directly to suppress import warnings
.runjags__findjags <- function() runjags::findjags()

# adapted from the runjags package version 2.2.0
RoBSA.private <- new.env()
# Use 'expression' for functions to avoid having to evaluate before the package is fully loaded:
assign("JAGS_path",       expression(.runjags__findjags()),           envir = RoBSA.private)
assign("RoBSA_version",   utils::packageVersion("RoBSA"),             envir = RoBSA.private)
assign("min_jags_major",  4,                                          envir = RoBSA.private)
assign("max_jags_major",  4,                                          envir = RoBSA.private)
assign("max_cores",       parallel::detectCores(logical = TRUE) - 1,  envir = RoBSA.private)
assign("check_scaling",   TRUE,                                       envir = RoBSA.private)
assign("distributions",   .distributions,                             envir = RoBSA.private)

# each worker fits whole models, one chain per core; machines with fewer cores
# than chains still get one worker
.parallel_workers <- function(max_cores, chains){
  return(max(1, floor(max_cores / chains)))
}

# check proper BayesTools package version
.check_BayesTools <- function(
    RoBSA.version      = try(utils::packageVersion("RoBSA")),
    BayesTools.version = try(utils::packageVersion("BayesTools"))){

  if(inherits(RoBSA.version, "try-error") | inherits(BayesTools.version, "try-error")){
    return(invisible(FALSE))
  }

  if(is.null(RoBSA.version) | is.null(BayesTools.version)){
    return(invisible(FALSE))
  }

  # a development version (x.y.z.9001) requires what its release requires
  RoBSA.release <- paste0(unlist(RoBSA.version)[1:3], collapse = ".")

  BayesTools_required <- switch(
    RoBSA.release,
    "1.0.0" = c("0.2.10", "0.2.13"),
    "1.0.1" = c("0.2.10", "0.2.13"),
    "1.0.2" = c("0.2.14", "0.2.23"),
    "1.0.3" = c("0.2.14", "0.2.23"),
    "1.0.4" = c("0.3.0", "999.999.999"),
    "1.1.0" = c("0.4.0", "999.999.999"),
    stop("New RoBSA version needs to be defined in '.check_BayesTools' function!")
  )

  min_OK <- BayesTools.version >= BayesTools_required[1]
  max_OK <- BayesTools.version <= BayesTools_required[2]

  if(min_OK && max_OK){
    return(invisible(TRUE))
  }else{
    warning(sprintf(
      "RoBSA version %1$s requires BayesTools version higher or equal %2$s and lower or equal %3$s.",
      paste0(RoBSA.version, collapse = "."),
      BayesTools_required[1],
      BayesTools_required[2]
    ), call.=FALSE)
    return(invisible(FALSE))
  }
}

