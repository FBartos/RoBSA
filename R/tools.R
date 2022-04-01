
### tools
.has_aux         <- function(distributions){
  distributions != "exp-aft"
}
.intercept_name  <- function(distribution){
  switch(
    distribution,
    "exp-aft"     = "rate",
    "weibull-aft" = "scale",
    "lnorm-aft"   = "meanlog",
    "llogis-aft"  = "scale",
    "gamma-aft"   = "rate"
  )
}
.aux_name        <- function(distribution){
  switch(
    distribution,
    "exp-aft"     = NULL,
    "weibull-aft" = "shape",
    "lnorm-aft"   = "sdlog",
    "llogis-aft"  = "shape",
    "gamma-aft"   = "shape"
  )
}
.intercept_transformation <- function(distribution){
  switch(
    distribution,
    "exp-aft"     = function(x)exp(-x),
    "weibull-aft" = function(x)exp(x),
    "lnorm-aft"   = function(x)x,
    "llogis-aft"  = function(x)exp(x),
    "gamma-aft"   = function(x)exp(-x)
  )
}
.aux_transformation       <- function(distribution){
  switch(
    distribution,
    "exp-aft"     = NULL,
    "weibull-aft" = function(x)x,
    "lnorm-aft"   = function(x)x,
    "llogis-aft"  = function(x)x,
    "gamma-aft"   = function(x)x
  )
}
.is_zero_spike   <- function(prior){
  if(prior$distribution == "point"){
    if(prior$parameters$location == 0){
      return(TRUE)
    }
  }
  return(FALSE)
}
.is_model_constant         <- function(priors){
  # checks whether there is at least one non-nill prior
  return(all(sapply(priors, function(prior) is.prior.point(prior) | is.prior.none(prior))))
}
