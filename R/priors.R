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
prior_factor <- function(distribution, parameters, truncation = list(lower = -Inf, upper = Inf), prior_weights = 1, contrast = "orthonormal"){
  BayesTools::prior_factor(distribution = distribution, parameters = parameters, truncation = truncation, prior_weights = prior_weights, contrast = contrast)
}

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


#' @title Mean difference contrast matrix
#'
#' @description Return a matrix of mean difference contrasts.
#' This is an adjustment to the \code{contr.orthonormal} that ascertains that the prior
#' distributions on difference between the gran mean and factor level are identical independent
#' of the number of factor levels (which does not hold for the orthonormal contrast). Furthermore,
#' the contrast is re-scaled so the specified prior distribution exactly corresponds to the prior
#' distribution on difference between each factor level and the grand mean -- this is approximately
#' twice the scale of \code{contr.orthonormal}.
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.meandif(c(1, 2))
#' contr.meandif(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
#' @export
contr.meandif <- BayesTools::contr.meandif
