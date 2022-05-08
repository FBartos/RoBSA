# some general functions
.distributions <- c(
  "exp-aft",
  "weibull-aft",
  "lnorm-aft",
  "llogis-aft",
  "gamma-aft"
)

.get_marginal_distribution_function <- function(distribution, type){
  eval(parse(text = paste0(".", gsub("-", "_", distribution))))[[paste0("log_", switch(
    type,
    "event"   = "density",
    "cens_r"  = "survival",
    "cens_l"  = stop("not implemented"),
    "cens_ir" = stop("not implemented"),
    "cens_il" = stop("not implemented")
  ))]]
}
.get_prediction_distribution_function <- function(distribution, type){
  eval(parse(text = paste0(".", gsub("-", "_", distribution))))[[type]]
}


#' @title Exponential AFT parametric family.
#'
#' @description (log) density, hazard, and survival
#' functions for AFT exponential parametric family.
#'
#' @param t vector of survival times
#' @param eta linear predictor
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#'
#'
#' @return \code{exp_aft_density}, \code{exp_aft_hazard}, and
#' \code{exp_aft_survival} return the density, hazard, and survival
#' of the specified survival distribution. The \code{exp_aft_log_density},
#' \code{exp_aft_log_hazard}, \code{exp_aft_log_survival} return log of
#' the corresponding qualities. \code{exp_aft_mean} and \code{exp_aft_sd}
#' return the mean and standard deviation of the specified survival distribution.
#' \code{exp_aft_r}, \code{exp_aft_q}, and \code{exp_aft_p} return a random
#' generation, quantiles, and cumulative probabilities of the specified
#' survival distribution.
#'
#'
#' @export exp_aft_log_density
#' @export exp_aft_log_hazard
#' @export exp_aft_log_survival
#' @export exp_aft_density
#' @export exp_aft_hazard
#' @export exp_aft_survival
#' @export exp_aft_mean
#' @export exp_aft_sd
#' @export exp_aft_r
#' @export exp_aft_q
#' @export exp_aft_p
#' @name exp-aft
NULL

# corresponds to stats::dexp(x, rate = 1/exp(eta))
.exp_aft <- list(
  log_density   = function(t, eta) .exp_aft$log_hazard(t, eta) + .exp_aft$log_survival(t, eta),
  log_hazard    = function(t, eta) -eta,
  log_survival  = function(t, eta) -t * exp(-eta),
  density       = function(t, eta) exp(.exp_aft$log_density(t, eta)),
  hazard        = function(t, eta) exp(.exp_aft$log_hazard(t, eta)),
  survival      = function(t, eta) exp(.exp_aft$log_survival(t, eta)),
  mean          = function(eta)    exp(eta),
  sd            = function(eta)    exp(eta),
  r             = function(n, eta) stats::rexp(n, rate = 1/exp(eta)),
  q             = function(p, eta) stats::qexp(p, rate = 1/exp(eta)),
  p             = function(q, eta) stats::pexp(q, rate = 1/exp(eta))
)

#' @rdname exp-aft
exp_aft_log_density  <- function(t, eta){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.exp_aft$log_density(t, eta))
}
#' @rdname exp-aft
exp_aft_log_hazard   <- function(t, eta){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.exp_aft$log_hazard(t, eta))
}
#' @rdname exp-aft
exp_aft_log_survival <- function(t, eta){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.exp_aft$log_survival(t, eta))
}
#' @rdname exp-aft
exp_aft_density      <- function(t, eta){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.exp_aft$density(t, eta))
}
#' @rdname exp-aft
exp_aft_hazard       <- function(t, eta){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.exp_aft$hazard(t, eta))
}
#' @rdname exp-aft
exp_aft_survival     <- function(t, eta){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.exp_aft$survival(t, eta))
}
#' @rdname exp-aft
exp_aft_mean <- function(eta){
  return(.exp_aft$mean(eta))
}
#' @rdname exp-aft
exp_aft_sd   <- function(eta){
  return(.exp_aft$sd(eta))
}
#' @rdname exp-aft
exp_aft_r    <- function(n, eta){
  BayesTools::check_int(n, "n", lower = 0)
  return(.exp_aft$r(n, eta))
}
#' @rdname exp-aft
exp_aft_q    <- function(p, eta){
  BayesTools::check_real(p, "p", lower = 0, upper = 1)
  return(.exp_aft$q(p, eta))
}
#' @rdname exp-aft
exp_aft_p    <- function(q, eta){
  BayesTools::check_real(q, "q", lower = 0)
  return(.exp_aft$p(q, eta))
}


#' @title Weibull AFT parametric family.
#'
#' @description (log) density, hazard, and survival
#' functions for AFT Weibull parametric family.
#'
#' @param t vector of survival times
#' @param eta linear predictor
#' @param shape auxiliary parameter
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#'
#'
#' @return \code{weibull_aft_density}, \code{weibull_aft_hazard}, and
#' \code{weibull_aft_survival} return the density, hazard, and survival
#' of the specified survival distribution. The \code{weibull_aft_log_density},
#' \code{weibull_aft_log_hazard}, \code{weibull_aft_log_survival} return log of
#' the corresponding qualities. \code{weibull_aft_mean} and \code{weibull_aft_sd}
#' return the mean and standard deviation of the specified survival distribution.
#' \code{weibull_aft_r}, \code{weibull_aft_q}, and \code{weibull_aft_p} return a random
#' generation, quantiles, and cumulative probabilities of the specified
#' survival distribution.
#'
#'
#' @export weibull_aft_log_density
#' @export weibull_aft_log_hazard
#' @export weibull_aft_log_survival
#' @export weibull_aft_density
#' @export weibull_aft_hazard
#' @export weibull_aft_survival
#' @export weibull_aft_mean
#' @export weibull_aft_sd
#' @export weibull_aft_r
#' @export weibull_aft_q
#' @export weibull_aft_p
#' @name weibull-aft
NULL

# corresponds to stats::dweibull(x, shape = shape, scale = exp(eta))
.weibull_aft <- list(
  log_density   = function(t, eta, shape) .weibull_aft$log_hazard(t, eta, shape) + .weibull_aft$log_survival(t, eta, shape),
  log_hazard    = function(t, eta, shape) log(shape) + (shape - 1) * log(t) - shape * eta,
  log_survival  = function(t, eta, shape) -t^shape * exp(-shape * eta),
  density       = function(t, eta, shape) exp(.weibull_aft$log_density(t, eta, shape)),
  hazard        = function(t, eta, shape) exp(.weibull_aft$log_hazard(t, eta, shape)),
  survival      = function(t, eta, shape) exp(.weibull_aft$log_survival(t, eta, shape)),
  mean          = function(eta, shape)    exp(eta)*gamma(1+1/shape),
  sd            = function(eta, shape)    sqrt(exp(eta)^2 * (gamma(1+2/shape) - gamma(1+1/shape)^2)),
  r             = function(n, eta, shape) stats::rweibull(n, shape = shape, scale = exp(eta)),
  q             = function(p, eta, shape) stats::qweibull(p, shape = shape, scale = exp(eta)),
  p             = function(q, eta, shape) stats::pweibull(q, shape = shape, scale = exp(eta))
)

#' @rdname weibull-aft
weibull_aft_log_density  <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.weibull_aft$log_density(t, eta, shape))
}
#' @rdname weibull-aft
weibull_aft_log_hazard   <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.weibull_aft$log_hazard(t, eta, shape))
}
#' @rdname weibull-aft
weibull_aft_log_survival <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.weibull_aft$log_survival(t, eta, shape))
}
#' @rdname weibull-aft
weibull_aft_density      <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.weibull_aft$density(t, eta, shape))
}
#' @rdname weibull-aft
weibull_aft_hazard       <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.weibull_aft$hazard(t, eta, shape))
}
#' @rdname weibull-aft
weibull_aft_survival     <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.weibull_aft$survival(t, eta, shape))
}
#' @rdname weibull-aft
weibull_aft_mean <- function(eta, shape){
  return(.weibull_aft$mean(eta, shape))
}
#' @rdname weibull-aft
weibull_aft_sd   <- function(eta, shape){
  return(.weibull_aft$sd(eta, shape))
}
#' @rdname weibull-aft
weibull_aft_r    <- function(n, eta, shape){
  BayesTools::check_int(n, "n", lower = 0)
  return(.weibull_aft$r(n, eta, shape))
}
#' @rdname weibull-aft
weibull_aft_q    <- function(p, eta, shape){
  BayesTools::check_real(p, "p", lower = 0, upper = 1)
  return(.weibull_aft$q(p, eta, shape))
}
#' @rdname weibull-aft
weibull_aft_p    <- function(q, eta, shape){
  BayesTools::check_real(q, "q", lower = 0)
  return(.weibull_aft$p(q, eta, shape))
}


#' @title Log-normal AFT parametric family.
#'
#' @description (log) density, hazard, and survival
#' functions for AFT log-normal parametric family.
#'
#' @param t vector of survival times
#' @param eta linear predictor
#' @param sd auxiliary parameter
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#'
#'
#' @return \code{lnorm_aft_density}, \code{lnorm_aft_hazard}, and
#' \code{lnorm_aft_survival} return the density, hazard, and survival
#' of the specified survival distribution. The \code{lnorm_aft_log_density},
#' \code{lnorm_aft_log_hazard}, \code{lnorm_aft_log_survival} return log of
#' the corresponding qualities. \code{lnorm_aft_mean} and \code{lnorm_aft_sd}
#' return the mean and standard deviation of the specified survival distribution.
#' \code{lnorm_aft_r}, \code{lnorm_aft_q}, and \code{lnorm_aft_p} return a random
#' generation, quantiles, and cumulative probabilities of the specified
#' survival distribution.
#'
#'
#' @export lnorm_aft_log_density
#' @export lnorm_aft_log_hazard
#' @export lnorm_aft_log_survival
#' @export lnorm_aft_density
#' @export lnorm_aft_hazard
#' @export lnorm_aft_survival
#' @export lnorm_aft_mean
#' @export lnorm_aft_sd
#' @export lnorm_aft_r
#' @export lnorm_aft_q
#' @export lnorm_aft_p
#' @name lnorm-aft
NULL

# corresponds to stats::dlnorm(x, meanlog = eta, sdlog = sd)
.lnorm_aft <- list(
  log_density   = function(t, eta, sd) stats::dlnorm(t, meanlog = eta, sdlog = sd, log = TRUE),
  log_hazard    = function(t, eta, sd) .lnorm_aft$log_density(t, eta, sd) - .lnorm_aft$log_survival(t, eta, sd),
  log_survival  = function(t, eta, sd) stats::plnorm(t, meanlog = eta, sdlog = sd, lower.tail = FALSE, log.p = TRUE),
  density       = function(t, eta, sd) exp(.lnorm_aft$log_density(t, eta, sd)),
  hazard        = function(t, eta, sd) exp(.lnorm_aft$log_hazard(t, eta, sd)),
  survival      = function(t, eta, sd) exp(.lnorm_aft$log_survival(t, eta, sd)),
  mean          = function(eta, sd)    exp(eta + sd^2/2),
  sd            = function(eta, sd)    sqrt( (exp(sd^2)-1) * exp(2*eta+sd^2) ),
  r             = function(n, eta, sd) stats::rlnorm(n, meanlog = eta, sdlog = sd),
  q             = function(p, eta, sd) stats::qlnorm(p, meanlog = eta, sdlog = sd),
  p             = function(q, eta, sd) stats::plnorm(q, meanlog = eta, sdlog = sd)
)

#' @rdname lnorm-aft
lnorm_aft_log_density  <- function(t, eta, sd){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.lnorm_aft$log_density(t, eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_log_hazard   <- function(t, eta, sd){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.lnorm_aft$log_hazard(t, eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_log_survival <- function(t, eta, sd){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.lnorm_aft$log_survival(t, eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_density      <- function(t, eta, sd){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.lnorm_aft$density(t, eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_hazard       <- function(t, eta, sd){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.lnorm_aft$hazard(t, eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_survival     <- function(t, eta, sd){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.lnorm_aft$survival(t, eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_mean <- function(eta, sd){
  return(.lnorm_aft$mean(eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_sd   <- function(eta, sd){
  return(.lnorm_aft$sd(eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_r    <- function(n, eta, sd){
  BayesTools::check_int(n, "n", lower = 0)
  return(.lnorm_aft$r(n, eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_q    <- function(p, eta, sd){
  BayesTools::check_real(p, "p", lower = 0, upper = 1)
  return(.lnorm_aft$q(p, eta, sd))
}
#' @rdname lnorm-aft
lnorm_aft_p    <- function(q, eta, sd){
  BayesTools::check_real(q, "q", lower = 0)
  return(.lnorm_aft$p(q, eta, sd))
}

#' @title Log-logistic AFT parametric family.
#'
#' @description (log) density, hazard, and survival
#' functions for AFT log-logistic parametric family.
#'
#' @param t vector of survival times
#' @param eta linear predictor
#' @param shape auxiliary parameter
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#'
#'
#' @return \code{llogis_aft_density}, \code{llogis_aft_hazard}, and
#' \code{llogis_aft_survival} return the density, hazard, and survival
#' of the specified survival distribution. The \code{llogis_aft_log_density},
#' \code{llogis_aft_log_hazard}, \code{llogis_aft_log_survival} return log of
#' the corresponding qualities. \code{llogis_aft_mean} and \code{llogis_aft_sd}
#' return the mean and standard deviation of the specified survival distribution.
#' \code{llogis_aft_r}, \code{llogis_aft_q}, and \code{llogis_aft_p} return a random
#' generation, quantiles, and cumulative probabilities of the specified
#' survival distribution.
#'
#'
#' @export llogis_aft_log_density
#' @export llogis_aft_log_hazard
#' @export llogis_aft_log_survival
#' @export llogis_aft_density
#' @export llogis_aft_hazard
#' @export llogis_aft_survival
#' @export llogis_aft_mean
#' @export llogis_aft_sd
#' @export llogis_aft_r
#' @export llogis_aft_q
#' @export llogis_aft_p
#' @name llogis-aft
NULL

# corresponds to flexsurv::dllogis(x, shape = shape, scale = exp(eta))
.llogis_aft <- list(
  log_density   = function(t, eta, shape) log(shape) - eta + (shape - 1) * (log(t) - eta) - 2 * log(1 + (t / exp(eta))^shape),
  log_hazard    = function(t, eta, shape) .llogis_aft$log_density(t, eta, shape) - .llogis_aft$log_survival(t, eta, shape),
  log_survival  = function(t, eta, shape) - log(1 + (t / exp(eta))^shape),
  density       = function(t, eta, shape) exp(.llogis_aft$log_density(t, eta, shape)),
  hazard        = function(t, eta, shape) exp(.llogis_aft$log_hazard(t, eta, shape)),
  survival      = function(t, eta, shape) exp(.llogis_aft$log_survival(t, eta, shape)),
  mean          = function(eta, shape)    ifelse(shape > 1, (exp(eta) * pi / shape) / sin (pi / shape), NA),
  sd            = function(eta, shape)    ifelse(shape > 2, sqrt(exp(eta)^2 * (2*pi/shape / sin(2*pi/shape) - (pi/shape)^2 / sin(pi/shape)^2)), NA),
  r             = function(n, eta, shape) flexsurv::rllogis(n, shape = shape, scale = exp(eta)),
  q             = function(p, eta, shape) flexsurv::qllogis(p, shape = shape, scale = exp(eta)),
  p             = function(q, eta, shape) flexsurv::pllogis(q, shape = shape, scale = exp(eta))
)

#' @rdname llogis-aft
llogis_aft_log_density  <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.llogis_aft$log_density(t, eta, shape))
}
#' @rdname llogis-aft
llogis_aft_log_hazard   <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.llogis_aft$log_hazard(t, eta, shape))
}
#' @rdname llogis-aft
llogis_aft_log_survival <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.llogis_aft$log_survival(t, eta, shape))
}
#' @rdname llogis-aft
llogis_aft_density      <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.llogis_aft$density(t, eta, shape))
}
#' @rdname llogis-aft
llogis_aft_hazard       <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.llogis_aft$hazard(t, eta, shape))
}
#' @rdname llogis-aft
llogis_aft_survival     <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.llogis_aft$survival(t, eta, shape))
}
#' @rdname llogis-aft
llogis_aft_mean <- function(eta, shape){
  if(any(shape <= 1))
    warning("llogis-aft: mean is undefinied for shape <= 1")
  return(.llogis_aft$mean(eta, shape))
}
#' @rdname llogis-aft
llogis_aft_sd   <- function(eta, shape){
  if(any(shape <= 2))
    warning("llogis-aft: sd is undefinied for shape <= 2")
  return(.llogis_aft$sd(eta, shape))
}
#' @rdname llogis-aft
llogis_aft_r    <- function(n, eta, shape){
  BayesTools::check_int(n, "n", lower = 0)
  return(.llogis_aft$r(n, eta, shape))
}
#' @rdname llogis-aft
llogis_aft_q    <- function(p, eta, shape){
  BayesTools::check_real(p, "p", lower = 0, upper = 1)
  return(.llogis_aft$q(p, eta, shape))
}
#' @rdname llogis-aft
llogis_aft_p    <- function(q, eta, shape){
  BayesTools::check_real(q, "q", lower = 0)
  return(.llogis_aft$p(q, eta, shape))
}


#' @title Gamma AFT parametric family.
#'
#' @description (log) density, hazard, and survival
#' functions for AFT gamma parametric family.
#'
#' @param t vector of survival times
#' @param eta linear predictor
#' @param shape auxiliary parameter
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#'
#'
#' @return \code{gamma_aft_density}, \code{gamma_aft_hazard}, and
#' \code{gamma_aft_survival} return the density, hazard, and survival
#' of the specified survival distribution. The \code{gamma_aft_log_density},
#' \code{gamma_aft_log_hazard}, \code{gamma_aft_log_survival} return log of
#' the corresponding qualities. \code{gamma_aft_mean} and \code{gamma_aft_sd}
#' return the mean and standard deviation of the specified survival distribution.
#' \code{gamma_aft_r}, \code{gamma_aft_q}, and \code{gamma_aft_p} return a random
#' generation, quantiles, and cumulative probabilities of the specified
#' survival distribution.
#'
#'
#' @export gamma_aft_log_density
#' @export gamma_aft_log_hazard
#' @export gamma_aft_log_survival
#' @export gamma_aft_density
#' @export gamma_aft_hazard
#' @export gamma_aft_survival
#' @export gamma_aft_mean
#' @export gamma_aft_sd
#' @export gamma_aft_r
#' @export gamma_aft_q
#' @export gamma_aft_p
#' @name gamma-aft
NULL

# corresponds to stats::dgamma(x, shape = shape, scale = 1/ exp(-eta))
.gamma_aft <- list(
  log_density   = function(t, eta, shape) stats::dgamma(t, shape = shape, scale = 1/ exp(-eta), log = TRUE),
  log_hazard    = function(t, eta, shape) .gamma_aft$log_density(t, eta, shape) - .gamma_aft$log_survival(t, eta, shape),
  log_survival  = function(t, eta, shape) stats::pgamma(t, shape = shape, scale = 1/ exp(-eta), lower.tail = FALSE, log.p = TRUE),
  density       = function(t, eta, shape) exp(.gamma_aft$log_density(t, eta, shape)),
  hazard        = function(t, eta, shape) exp(.gamma_aft$log_hazard(t, eta, shape)),
  survival      = function(t, eta, shape) exp(.gamma_aft$log_survival(t, eta, shape)),
  mean          = function(eta, shape)    shape / exp(-eta),
  sd            = function(eta, shape)    sqrt(shape / exp(-eta)^2),
  r             = function(n, eta, shape) stats::rgamma(n, shape = shape, scale = 1/exp(-eta)),
  q             = function(p, eta, shape) stats::qgamma(p, shape = shape, scale = 1/exp(-eta)),
  p             = function(q, eta, shape) stats::pgamma(q, shape = shape, scale = 1/exp(-eta))
)

#' @rdname gamma-aft
gamma_aft_log_density  <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.gamma_aft$log_density(t, eta, shape))
}
#' @rdname gamma-aft
gamma_aft_log_hazard   <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.gamma_aft$log_hazard(t, eta, shape))
}
#' @rdname gamma-aft
gamma_aft_log_survival <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.gamma_aft$log_survival(t, eta, shape))
}
#' @rdname gamma-aft
gamma_aft_density      <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.gamma_aft$density(t, eta, shape))
}
#' @rdname gamma-aft
gamma_aft_hazard       <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.gamma_aft$hazard(t, eta, shape))
}
#' @rdname gamma-aft
gamma_aft_survival     <- function(t, eta, shape){
  BayesTools::check_real(t, "t", lower = 0, allow_bound = TRUE, check_length = FALSE)
  return(.gamma_aft$survival(t, eta, shape))
}
#' @rdname gamma-aft
gamma_aft_mean <- function(eta, shape){
  return(.gamma_aft$mean(eta, shape))
}
#' @rdname gamma-aft
gamma_aft_sd   <- function(eta, shape){
  return(.gamma_aft$sd(eta, shape))
}
#' @rdname gamma-aft
gamma_aft_r    <- function(n, eta, shape){
  BayesTools::check_int(n, "n", lower = 0)
  return(.gamma_aft$r(n, eta, shape))
}
#' @rdname gamma-aft
gamma_aft_q    <- function(p, eta, shape){
  BayesTools::check_real(p, "p", lower = 0, upper = 1)
  return(.gamma_aft$q(p, eta, shape))
}
#' @rdname gamma-aft
gamma_aft_p    <- function(q, eta, shape){
  BayesTools::check_real(q, "q", lower = 0)
  return(.gamma_aft$p(q, eta, shape))
}


# some helper functions ----
lnorm_mean     <- function(meanlog, sdlog){
  exp(meanlog + sdlog^2/2)
}
lnorm_sd       <- function(meanlog, sdlog){
  sqrt( (exp(sdlog^2)-1) * exp(2*meanlog+sdlog^2) )
}
lnorm_meanlog  <- function(mean, sd){
  log(mean) - log(-(-mean^2-sd^2)/mean^2)/2
}
lnorm_sdlog    <- function(mean, sd){
  sqrt(log(-(-mean^2-sd^2)/mean^2))
}
dlnorm_mean_sd <- function(x, mean, sd){
  stats::dlnorm(x, lnorm_meanlog(mean, sd), lnorm_sdlog(mean, sd))
}
rlnorm_mean_sd <- function(n, mean, sd){
  stats::rlnorm(n, lnorm_meanlog(mean, sd), lnorm_sdlog(mean, sd))
}



