#' RoBSA: Robust Bayesian survival analysis
#'
#' A framework for estimating ensembles of parametric survival models
#' and using Bayesian model averaging to combine them. The
#' ensembles use Bayes factors to test for the presence or absence of the
#' individual components (e.g., treatment effect, parametric family, etc...) and
#' model-averages parameter estimates based on posterior model probabilities.
#' The user can define a wide range of informative priors for all parameters
#' of interest. The package provides convenient functions for summary, visualizations,
#' fit diagnostics, and prior distribution calibration.
#'
#' \code{\link{RoBSA}} fits parametric models ....
#'
#'
#' @name RoBSA-package
#' @author František Bartoš \email{f.bartos96@@gmail.com}
#' @keywords package
#' @aliases RoBSA-package RoBSA_package RoBSA.package
#' @docType package
#' @section User guide: The \bold{RoBSA user guide} vignette explains the
#' methods in detail, and gives several worked examples.  A further vignette
#' \bold{RoBSA-examples} gives a few more complicated examples, and users
#' are encouraged to submit their own.
#' @references
#' @importFrom survival Surv
"_PACKAGE"


#' @name Surv
#' @inherit survival::Surv
#' @export
Surv <- survival::Surv
