# Capture the build configuration at installation for the runtime ABI check.
.RoBSA_jags_build <- rjags::pkg.compile.flags("ALL")

.onLoad <- function(libname, pkgname){

  RoBSA.private$RoBSA_version <- utils::packageDescription(pkgname, fields = 'Version')
  RoBSA.private$lib_name <- libname
  .load_RoBSA_module(pkgname)

  setopts <- mget('.RoBSA.options', envir=.GlobalEnv, ifnotfound = list(.RoBSA.options = NULL))[[1]]
  if(!is.null(setopts)){
    if(!is.list(setopts)){
      warning('Ignoring invalid (non-list) specification for .RoBSA.options on loading the RoBSA package', call.=FALSE)
    }else{
      do.call('RoBSA.options', args = setopts)
    }
  }

  .check_BayesTools()

}

.onUnload <- function(libpath){
  # Finalize unused models while their distributions are still registered.
  gc()
  if("RoBSA" %in% rjags::list.modules()){
    rjags::unload.module("RoBSA", quiet = TRUE)
  }
}

.load_RoBSA_module <- function(pkgname = "RoBSA"){

  libname <- RoBSA.private$lib_name
  dll <- paste0(pkgname, .Platform$dynlib.ext)
  module_location <- file.path(libname, pkgname, "libs", .Platform$r_arch)
  if(file.exists(file.path(module_location, dll))){
    rjags::load.pkg.module(pkgname, libname, .RoBSA_jags_build, quiet = TRUE)
  }else{
    # pkgload compiles the development module in the source tree.
    src <- file.path(getNamespaceInfo(pkgname, "path"), "src")
    paths <- c(file.path(src, .Platform$r_arch), src)
    paths <- paths[file.exists(file.path(paths, dll))]
    if(!length(paths)){
      stop("The RoBSA JAGS module is missing. Reinstall RoBSA using JAGS 5 or later.", call. = FALSE)
    }
    module_location <- paths[[1L]]
    rjags::load.pkg.module(pkgname, libname, .RoBSA_jags_build, check.only = TRUE)
    rjags::load.module(pkgname, path = module_location, quiet = TRUE)
  }
  RoBSA.private$module_location <- normalizePath(module_location, winslash = "/")
  invisible(TRUE)
}
