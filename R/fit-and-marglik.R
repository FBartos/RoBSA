




.fitting_priority  <- function(models){

  terms         <- sapply(models, function(m)sum(!sapply(m[["priors"]][["terms"]], BayesTools::is.prior.point)))
  distributions <- sapply(models, function(m)switch(
    as.character(m$distribution),
    "exp-aft"     = 0,
    "weibull-aft" = 1.5,
    "lnorm-aft"   = 1,
    "llogis-aft"  = 1,
    "gamma-aft"   = 1.5
  ))

  fitting_difficulty <- terms + distributions

  return(order(fitting_difficulty, decreasing = TRUE))
}
