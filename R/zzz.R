
.onLoad <- function(libname, pkgname){

  # load runjags
  requireNamespace("runjags")
  requireNamespace("survival")

  hereIsTheModule <- file.path(libname, pkgname)
  path <- file.path(hereIsTheModule, paste0("libs", Sys.getenv("R_ARCH")))
  tryCatch(rjags::load.module("RoBSA", path = path), error = function(e) warning(sprintf("The RoBSA module couldn't be loaded from %s. libname: %s, pkgname: %s.\n", path, libname, pkgname)))

}

.onAttach <- function(libname, pkgname){

  packageStartupMessage(
    "Hello :)")

}
