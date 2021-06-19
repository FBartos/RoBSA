# from RoBMA
.marglik_fail          <- function(){
  marg_lik        <- NULL
  marg_lik$logml  <- -Inf
  class(marg_lik) <- "bridge"
  return(marg_lik)
}
.marglik_distribution  <- function(x, par_name, distribution, data, truncation){

  if(distribution == "normal"){
    log_lik <- stats::dnorm(x, mean = data[[paste0("prior_",par_name,"_mean")]], sd = data[[paste0("prior_",par_name,"_sd")]], log = TRUE) -
      log(
        stats::pnorm(truncation$upper, data[[paste0("prior_",par_name,"_mean")]], data[[paste0("prior_",par_name,"_sd")]], lower.tail = TRUE, log.p = FALSE) -
          stats::pnorm(truncation$lower, data[[paste0("prior_",par_name,"_mean")]], data[[paste0("prior_",par_name,"_sd")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "lognormal"){
    log_lik <- stats::dlnorm(x, meanlog = data[[paste0("prior_",par_name,"_meanlog")]], sdlog = data[[paste0("prior_",par_name,"_sdlog")]], log = TRUE) -
      log(
        stats::plnorm(truncation$upper, data[[paste0("prior_",par_name,"_meanlog")]], data[[paste0("prior_",par_name,"_sdlog")]], lower.tail = TRUE, log.p = FALSE) -
          stats::plnorm(truncation$lower, data[[paste0("prior_",par_name,"_meanlog")]], data[[paste0("prior_",par_name,"_sdlog")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "t"){
    log_lik <- extraDistr::dlst(x, df = data[[paste0("prior_",par_name,"_df")]], mu = data[[paste0("prior_",par_name,"_location")]], sigma = data[[paste0("prior_",par_name,"_scale")]], log = TRUE) -
      log(
        extraDistr::plst(truncation$upper, df = data[[paste0("prior_",par_name,"_df")]], mu = data[[paste0("prior_",par_name,"_location")]], sigma = data[[paste0("prior_",par_name,"_scale")]], lower.tail = TRUE, log.p = FALSE) -
          extraDistr::plst(truncation$lower, df = data[[paste0("prior_",par_name,"_df")]], mu = data[[paste0("prior_",par_name,"_location")]], sigma = data[[paste0("prior_",par_name,"_scale")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "gamma"){
    log_lik <- stats::dgamma(x, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_rate")]], log = TRUE)  -
      log(
        stats::pgamma(truncation$upper, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_rate")]], lower.tail = TRUE, log.p = FALSE) -
          stats::pgamma(truncation$lower, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_rate")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "invgamma"){
    log_lik <- stats::dgamma(1/x, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_scale")]], log = TRUE) -
      log(
        stats::pgamma(truncation$lower^-1, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_scale")]], lower.tail = TRUE, log.p = FALSE) -
          stats::pgamma(truncation$upper^-1, shape = data[[paste0("prior_",par_name,"_shape")]], rate = data[[paste0("prior_",par_name,"_scale")]], lower.tail = TRUE, log.p = FALSE)
      )
  }else if(distribution == "uniform"){
    log_lik <- stats::dunif(x, min = data[[paste0("prior_",par_name,"_a")]], max = data[[paste0("prior_",par_name,"_b")]], log = TRUE)
  }else if(distribution == "point"){
    log_lik <- 0
  }

  return(log_lik)
}

