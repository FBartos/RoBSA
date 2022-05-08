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
#' @details Further details can be found in \insertCite{erp2017estimates;textual}{BayesTools},
#' \insertCite{gronau2017bayesian;textual}{BayesTools}, and
#' \insertCite{bartos2021bayesian;textual}{BayesTools}.
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

#' @title Orthornomal contrast matrix
#'
#' @description Return a matrix of orthornomal contrasts.
#' Code is based on \code{stanova::contr.bayes} and corresponding to description
#' by \insertCite{rouder2012default;textual}{BayesTools}
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.orthonormal(c(1, 2))
#' contr.orthonormal(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
#' @export
contr.orthonormal <- BayesTools::contr.orthonormal
