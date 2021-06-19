# corresponds to stats::dexp(x, rate = 1/exp(eta))
exp_aft_log_density  <- function(t, eta){
  exp_aft_log_hazard(t, eta) + exp_aft_log_survival(t, eta)
}
exp_aft_log_hazard   <- function(t, eta){
  -eta
}
exp_aft_log_survival <- function(t, eta){
  -t * exp(-eta)
}
exp_aft_density      <- function(t, eta){
  exp(exp_aft_log_density(t, eta))
}
exp_aft_hazard       <- function(t, eta){
  exp(exp_aft_log_hazard(t, eta))
}
exp_aft_survival     <- function(t, eta){
  exp(exp_aft_log_survival(t, eta))
}

exp_aft_mean <- function(eta){
  exp(eta)
}
exp_aft_sd   <- function(eta){
  exp(eta)
}
exp_aft_r    <- function(n, eta){
  rexp(n, rate = 1/exp(eta))
}
exp_aft_q    <- function(p, eta){
  qexp(p, rate = 1/exp(eta))
}
exp_aft_p    <- function(q, eta){
  pexp(q, rate = 1/exp(eta))
}


# corresponds to stats::dweibull(x, shape = shape, scale = exp(eta))
weibull_aft_log_density  <- function(t, eta, shape){
  weibull_aft_log_hazard(t, eta, shape) + weibull_aft_log_survival(t, eta, shape)
}
weibull_aft_log_hazard   <- function(t, eta, shape){
  log(shape) + (shape - 1) * log(t) - shape * eta
}
weibull_aft_log_survival <- function(t, eta, shape){
  -t^shape * exp(-shape * eta)
}
weibull_aft_density      <- function(t, eta, shape){
  exp(weibull_aft_log_density(t, eta, shape))
}
weibull_aft_hazard       <- function(t, eta, shape){
  exp(weibull_aft_log_hazard(t, eta, shape))
}
weibull_aft_survival     <- function(t, eta, shape){
  exp(weibull_aft_log_survival(t, eta, shape))
}

weibull_aft_mean <- function(eta, shape){
  exp(eta)*gamma(1+1/shape)
}
weibull_aft_sd   <- function(eta, shape){
  sqrt(exp(eta)^2 * (gamma(1+2/shape) - gamma(1+1/shape)^2))
}
weibull_aft_r    <- function(n, eta, shape){
  rweibull(n, shape = shape, scale = exp(eta))
}
weibull_aft_q    <- function(p, eta, shape){
  qweibull(p, shape = shape, scale = exp(eta))
}
weibull_aft_p    <- function(q, eta, shape){
  pweibull(q, shape = shape, scale = exp(eta))
}


# corresponds to stats::dlnorm(x, meanlog = eta, sdlog = sd)
lnorm_aft_log_density  <- function(t, eta, sd){
  dlnorm(t, meanlog = eta, sdlog = sd, log = TRUE)
}
lnorm_aft_log_hazard   <- function(t, eta, sd){
  lnorm_aft_log_density(t, eta, sd) - lnorm_aft_log_survival(t, eta, sd)
}
lnorm_aft_log_survival <- function(t, eta, sd){
  plnorm(t, meanlog = eta, sdlog = sd, lower.tail = FALSE, log.p = TRUE)
}
lnorm_aft_density      <- function(t, eta, sd){
  exp(lnorm_aft_log_density(t, eta, sd))
}
lnorm_aft_hazard       <- function(t, eta, sd){
  exp(lnorm_aft_log_hazard(t, eta, sd))
}
lnorm_aft_survival     <- function(t, eta, sd){
  exp(lnorm_aft_log_survival(t, eta, sd))
}

lnorm_aft_mean    <- function(eta, sd){
  exp(eta + sd^2/2)
}
lnorm_aft_sd      <- function(eta, sd){
  sqrt( (exp(sd^2)-1) * exp(2*eta+sd^2) )
}
lnorm_aft_r       <- function(n, eta, sd){
  rlnorm(n, meanlog = eta, sdlog = sd)
}
lnorm_aft_q       <- function(p, eta, sd){
  qlnorm(p, meanlog = eta, sdlog = sd)
}
lnorm_aft_p       <- function(q, eta, sd){
  plnorm(q, meanlog = eta, sdlog = sd)
}


# corresponds to flexsurv::dllogis(x, shape = shape, scale = exp(eta))
llogis_aft_log_density  <- function(t, eta, shape){
  log(shape) - eta + (shape - 1) * (log(t) - eta) - 2 * log(1 + (t / exp(eta))^shape)
}
llogis_aft_log_hazard   <- function(t, eta, shape){
  llogis_aft_log_density(t, eta, shape) - llogis_aft_log_survival(t, eta, shape)
}
llogis_aft_log_survival <- function(t, eta, shape){
  - log(1 + (t / exp(eta))^shape)
}
llogis_aft_density      <- function(t, eta, shape){
  exp(llogis_aft_log_density(t, eta, shape))
}
llogis_aft_hazard       <- function(t, eta, shape){
  exp(llogis_aft_log_hazard(t, eta, shape))
}
llogis_aft_survival     <- function(t, eta, shape){
  exp(llogis_aft_log_survival(t, eta, shape))
}

# only for shape > 1 and shape > 2
llogis_aft_mean <- function(eta, shape){
  if(any(shape <= 1))
    warning("mean is undefinied for shape <= 1")
  ifelse(shape > 1, (exp(eta) * pi / shape) / sin (pi / shape), NA)
}
llogis_aft_sd   <- function(eta, shape){
  if(any(shape <= 2))
    warning("sd is undefinied for shape <= 2")
  ifelse(shape > 2, sqrt(exp(eta)^2 * (2*pi/shape / sin(2*pi/shape) - (pi/shape)^2 / sin(pi/shape)^2)), NA)
}
llogis_aft_r    <- function(n, eta, shape){
  flexsurv::rllogis(n, shape = shape, scale = exp(eta))
}
llogis_aft_q    <- function(p, eta, shape){
  flexsurv::qllogis(p, shape = shape, scale = exp(eta))
}
llogis_aft_p    <- function(q, eta, shape){
  flexsurv::pllogis(q, shape = shape, scale = exp(eta))
}


# corresponds to stats::dgamma(x, shape = shape, scale = 1/ exp(-eta))
gamma_aft_log_density  <- function(t, eta, shape){
  dgamma(t, shape = shape, scale = 1/ exp(-eta), log = TRUE)
}
gamma_aft_log_hazard   <- function(t, eta, shape){
  gamma_aft_log_density(t, eta, shape) - gamma_aft_log_survival(t, eta, shape)
}
gamma_aft_log_survival <- function(t, eta, shape){
  pgamma(t, shape = shape, scale = 1/ exp(-eta), lower.tail = FALSE, log.p = TRUE)
}
gamma_aft_density      <- function(t, eta, shape){
  exp(gamma_aft_log_density(t, eta, shape))
}
gamma_aft_hazard       <- function(t, eta, shape){
  exp(gamma_aft_log_hazard(t, eta, shape))
}
gamma_aft_survival     <- function(t, eta, shape){
  exp(gamma_aft_log_survival(t, eta, shape))
}

gamma_aft_mean <- function(eta, shape){
  shape / exp(-eta)
}
gamma_aft_sd   <- function(eta, shape){
  sqrt(shape / exp(-eta)^2)
}
gamma_aft_r    <- function(n, eta, shape){
  rgamma(n, shape = shape, scale = 1/exp(-eta))
}
gamma_aft_q    <- function(p, eta, shape){
  qgamma(p, shape = shape, scale = 1/exp(-eta))
}
gamma_aft_p    <- function(q, eta, shape){
  pgamma(q, shape = shape, scale = 1/exp(-eta))
}


# some helper functions
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
  dlnorm(x, lnorm_meanlog(mean, sd), lnorm_sdlog(mean, sd))
}
rlnorm_mean_sd <- function(n, mean, sd){
  rlnorm(n, lnorm_meanlog(mean, sd), lnorm_sdlog(mean, sd))
}
