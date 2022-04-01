# adapted from the runjags package version 2.2.0
.onLoad <- function(libname, pkgname){

  requireNamespace("runjags")

  RoBSA.private$RoBSA_version <- utils::packageDescription(pkgname, fields = 'Version')

  # Get and save the library location, getting rid of any trailing / caused by r_arch being empty:
  module_location <- gsub('/$','', file.path(libname, pkgname, 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
  if(!file.exists(file.path(module_location, paste('RoBSA', .Platform$dynlib.ext, sep='')))){
    module_location <- NULL
    warning('The RoBSA module could not be loaded.', call. = FALSE)
  }else{
    rjags::load.module("RoBSA", path = module_location)
  }

  RoBSA.private$module_location <- module_location
  RoBSA.private$lib_name        <- libname

  setopts <- mget('.RoBSA.options', envir=.GlobalEnv, ifnotfound = list(.RoBSA.options = NULL))[[1]]
  if(!is.null(setopts)){
    if(!is.list(setopts)){
      warning('Ignoring invalid (non-list) specification for .RoBSA.options on loading the RoBSA package', call.=FALSE)
    }else{
      newopts <- do.call('RoBSA.options', args = setopts)
    }
  }

  .check_BayesTools()

}

.onAttach <- function(libname, pkgname){

  # nothing to do

}

.onUnload <- function(libpath){

  # tricking the dyn.library unload
  if(!is.null(RoBSA.private$lib_name)){
    library.dynam("RoBSA", "RoBSA", RoBSA.private$lib_name)
  }

  # Just in case it is not always safe to try and access an element of an env that is in the process of being deleted (when R quits):
  if(!is.null(RoBSA.private$module_location)){
    rjags::unload.module("RoBSA")
  }
}

.load_RoBSA_module <- function(pkgname = "RoBSA"){

  if(is.null(RoBSA.private$module_location) || (!is.null(RoBSA.private$module_location) && RoBSA.private$module_location == "")){
    libnames         <- .libPaths()
    module_locations <- sapply(libnames, function(libname) gsub('/$','', file.path(libname, pkgname, 'libs', if(.Platform$r_arch!="") .Platform$r_arch else "")))
    sapply(module_locations, function(module_location) rjags::load.module("RoBSA", path = module_location))
  }else{
    rjags::load.module("RoBSA", path = RoBSA.private$module_location)
  }

}
