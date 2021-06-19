##' RoBSA: Robust Bayesian survival analysis
##'
##' RoBSA: Bayesian model-averaged parametric survival models with informed prior
##' prior distributions and inclusion Bayes factors.
##'
##' \code{\link{RoBSA}} fits parametric models ....
##'
##'
##' @name RoBSA-package
##' @aliases RoBSA-package RoBSA
##' @docType package
##' @section User guide: The \bold{RoBSA user guide} vignette explains the
##' methods in detail, and gives several worked examples.  A further vignette
##' \bold{RoBSA-examples} gives a few more complicated examples, and users
##' are encouraged to submit their own.
##' @author František Bartoš \email{f.bartos96@@gmail.com}
##' @references
##' @keywords package
##' @useDynLib RoBSA, .registration = TRUE
##' @importFrom survival Surv

"_PACKAGE"

.onUnload <- function(libpath) {
  library.dynam.unload("RoBSA", libpath)
}
