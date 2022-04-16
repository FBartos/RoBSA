#' @name prior
#' @inherit BayesTools::prior
#' @export
prior <- BayesTools::prior

#' @name prior_none
#' @inherit BayesTools::prior_none
#' @export
prior_none <- BayesTools::prior_none

#' @name prior_factor
#' @inherit BayesTools::prior_factor
#' @export
prior_factor <- BayesTools::prior_factor

#' @name prior_informed
#' @inherit BayesTools::prior_informed
#' @details Further details can be found in \insertCite{erp2017estimates;textual}{RoBSA},
#' \insertCite{gronau2017bayesian;textual}{RoBSA}, and
#' \insertCite{bartos2021bayesian;textual}{RoBSA}.
#' @references
#' \insertAllCited{}
#' @export
prior_informed <- BayesTools::prior_informed

#' @name prior_informed_medicine_names
#' @title Names of medical subfields from the Cochrane database of systematic reviews
#' @description Contain names identifying the individual subfields
#' from the Cochrane database of systematic reviews. The individual
#' elements correspond to valid name arguments for the [prior_informed()]
#' function.
#' @export
prior_informed_medicine_names <- BayesTools::prior_informed_medicine_names

#' @name contr.orthonormal
#' @inherit BayesTools::contr.orthonormal
#' @export
contr.orthonormal <- BayesTools::contr.orthonormal
