#' @importFrom rjags load.pkg.module unload.module list.modules

# Package environment variables:
.jags_pkg_env <- new.env(parent = emptyenv())
.jags_pkg_env$pkgname <- character(0)
.jags_pkg_env$libname <- character(0)
.jags_pkg_env$dyn_loaded <- FALSE
.jags_pkg_env$pkg_compile_flags <- rjags::pkg.compile.flags("ALL")

# adapted from the runjags package version 2.2.0
.onLoad <- function(libname, pkgname){

  # Preserve the name of the package for use in .onUnload etc:
  .jags_pkg_env$pkgname <- pkgname
  .jags_pkg_env$libname <- libname

  # First try to load rjags (this sets paths on Windows etc):
  tryCatch({
    if(!isNamespaceLoaded("rjags")) loadNamespace("rjags")
  }, error = function(e) {
    stop("The following error was obtained when loading rjags:\n\t",
         conditionMessage(e), call.=FALSE)
  })
  # Ensure the version of rjags is compatible with what we need:
  rjags_version <- utils::packageVersion("rjags")
  if(rjags_version >= "5" && rjags_version < "5.4"){
    stop("An incompatible version of rjags (", gsub(".","-",rjags_version,fixed=TRUE), ") is installed\n\tPlease install either rjags 4-18 or rjags >= 5-4", call.=FALSE)
  }

  # Guard against non-present dynlib when running roxygen:
  if(file.exists(
    if(nzchar(.Platform$r_arch)){
      file.path(libname, pkgname, "libs", .Platform$r_arch, paste0(pkgname, .Platform$dynlib.ext))
    }else{
      file.path(libname, pkgname, "libs", paste0(pkgname, .Platform$dynlib.ext))
    }
  )){
    tryCatch({
      # First load the package module, so that rjags can report on any JAGS
      # version incompatibility:
      rjags::load.pkg.module(pkgname, libname, .jags_pkg_env$pkg_compile_flags, quiet=TRUE)

      # Then load the dynamic library:
      library.dynam(pkgname, pkgname, libname)
      # And set the dyn_loaded flag to TRUE:
      .jags_pkg_env$dyn_loaded <- TRUE

    }, error = function(e) {
      if(.jags_pkg_env$pkgname %in% list.modules()){
        unload.module(.jags_pkg_env$pkgname, quiet=TRUE)
      }
      stop("Loading the package module/dynlib failed with the following error:\n\t", conditionMessage(e), call.=FALSE)
    })
  }

  requireNamespace("runjags")

  RoBSA.private$RoBSA_version <- utils::packageDescription(pkgname, fields = 'Version')

  # Get and save the library location, getting rid of any trailing / caused by r_arch being empty:
  module_location <- gsub('/$','', file.path(libname, pkgname, 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
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

  # Ignore if the dynlib failed to load at startup:
  if(.jags_pkg_env$dyn_loaded){
    # First unload the JAGS module (silently), if loaded:
    if(.jags_pkg_env$pkgname %in% list.modules()){
      unload.module(.jags_pkg_env$pkgname, quiet = TRUE)
    }
    # Then the dynlib:
    library.dynam.unload(.jags_pkg_env$pkgname, libpath)
  }

}

.load_RoBSA_module <- function(pkgname = "RoBSA"){

  rjags::load.pkg.module(.jags_pkg_env$pkgname, .jags_pkg_env$libname, .jags_pkg_env$pkg_compile_flags, quiet=FALSE)

}
