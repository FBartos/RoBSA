# from RoBMA (with removal of some arguments)

.fit_inits_mu_tau          <- function(prior, par){

  temp_x  <- NULL

  # generate the value
  if(prior$distribution == "normal"){
    while(length(temp_x) != 1){
      temp_x <- stats::rnorm(1, mean = prior$parameters$mean, sd = prior$parameters$sd)
      temp_x <- temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper]
    }
  }else if(prior$distribution == "lognormal"){
    while(length(temp_x) != 1){
      temp_x <- stats::rlnorm(1, meanlog = prior$parameters$meanlog, sdlog = prior$parameters$sdlog)
      temp_x <- temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper]
    }
  }else if(prior$distribution == "t"){
    while(length(temp_x) != 1){
      temp_x <- extraDistr::rlst(1, df = prior$parameters$df, mu = prior$parameters$location, sigma = prior$parameters$scale)
      temp_x <- temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper]
    }
  }else if(prior$distribution == "gamma"){
    while(length(temp_x) != 1){
      temp_x <- stats::rgamma(1, shape = prior$parameters$shape, rate = prior$parameters$rate)
      temp_x <- temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper]
    }
  }else if(prior$distribution == "invgamma"){
    while(length(temp_x) != 1){
      temp_x <- stats::rgamma(1, shape = prior$parameters$shape, rate = prior$parameters$scale)
      temp_x <- temp_x[temp_x >= prior$truncation$upper^-1 & temp_x <= prior$truncation$lower^-1]
    }
  }else if(prior$distribution == "uniform"){
    temp_x <- stats::runif(1, min = prior$parameters$a, max = prior$parameters$b)
  }

  # name the parameter
  if(par %in% c("PET", "PEESE")){
    if(prior$distribution == "invgamma"){
      names(temp_x) <- paste0("inv_", par)
    }else if(prior$distribution %in% c("normal", "lognormal", "t", "gamma", "uniform")){
      names(temp_x) <- par
    }
  }else{
    if(prior$distribution == "invgamma"){
      names(temp_x) <- paste0("inv_", par)
    }else if(prior$distribution %in% c("normal", "lognormal", "t", "gamma", "uniform")){
      names(temp_x) <- par
    }
  }

  return(temp_x)
}
.JAGS_distribution         <- function(par_name, distribution, truncation, priors_scale, effect_measure){

  if(par_name %in% c("PET", "PEESE")){
    par_data <- "omega"
  }else{
    par_data <- par_name
  }

  par_name_transform <- paste0(par_name, "_transform")

  # distribution
  if(distribution == "point"){
    syntax <- paste0(par_name, " = prior_",par_data,"_location\n")
  }else if(distribution == "normal"){
    syntax <- paste0(par_name," ~ dnorm(prior_",par_data,"_mean, pow(prior_",par_data,"_sd, -2))")
  }else if(distribution == "lognormal"){
    syntax <- paste0(par_name," ~ dlnorm(prior_",par_data,"_meanlog, pow(prior_",par_data,"_sdlog, -2))")
  }else if(distribution == "t"){
    syntax <- paste0(par_name," ~ dt(prior_",par_data,"_location, pow(prior_",par_data,"_scale, -2), prior_",par_data,"_df)")
  }else if(distribution == "gamma"){
    syntax <- paste0(par_name," ~ dgamma(prior_",par_data,"_shape, prior_",par_data,"_rate)")
  }else if(distribution == "invgamma"){
    syntax <- paste0("inv_",par_name," ~ dgamma(prior_",par_data,"_shape, prior_",par_data,"_scale)")
  }else if(distribution == "uniform"){
    syntax <- paste0(par_name," ~ dunif(prior_",par_data,"_a, prior_",par_data,"_b)\n")
  }

  # truncation
  if(!distribution %in% c("point", "uniform")){
    if(!(is.infinite(truncation$lower)  & is.infinite(truncation$lower))){
      # the truncation for invgamma needs to be done the other way around since we sample from gamma
      if(distribution == "invgamma"){
        syntax <- paste0(syntax, "T(",
                         ifelse(is.infinite(truncation$upper^-1),"",truncation$upper^-1),
                         ",",
                         ifelse(is.infinite(truncation$lower^-1),"",truncation$lower^-1),
                         ")\n")
      }else{
        syntax <- paste0(syntax, "T(",
                         ifelse(is.infinite(truncation$lower),"",truncation$lower),
                         ",",
                         ifelse(is.infinite(truncation$upper),"",truncation$upper),
                         ")\n")
      }
    }else{
      syntax <- paste0(syntax, "\n")
    }
  }

  # transformations
  if(distribution == "invgamma"){
    syntax <- paste0(syntax, par_name," = pow(inv_",par_name,", -1)\n")
  }

  if(par_name %in% c("mu", "tau")){
    syntax <- paste0(syntax, .JAGS_scale(priors_scale, effect_measure, par_name, par_name_transform))
  }else if(par_name %in% c("PET")){
    syntax <- paste0(syntax, paste0(par_name_transform, " = ", par_name, "\n"))
  }else if(par_name %in% c("PEESE")){
    # don't forget that the transformation is inverse for PEESE
    syntax <- paste0(syntax, .JAGS_scale(effect_measure, priors_scale, par_name, par_name_transform))
  }

  # # the precise transformation are not used due the inability to rescale large variances
  # # instead, aproximate linear scaling is employed in the same way as in metaBMA package
  # if(par_name == "mu"){
  #   syntax <- paste0(syntax, .JAGS_transformation(priors_scale, effect_measure, par_name, par_name_transform))
  # }else if(par_name == "tau"){
  #   syntax <- paste0(syntax, .JAGS_transformation_se(priors_scale, effect_measure, par_name, "mu", par_name_transform))
  # }

  return(syntax)
}
.JAGS_distribution_omega   <- function(distribution, truncation, parameters, priors_scale, effect_measure){

  if(distribution == "point"){
    return()
  }else if(grepl("PET", distribution)){
    syntax <- .JAGS_distribution("PET", substr(distribution, 5, nchar(distribution)), truncation, priors_scale, effect_measure)
  }else if(grepl("PEESE", distribution)){
    syntax <- .JAGS_distribution("PEESE", substr(distribution, 7, nchar(distribution)), truncation, priors_scale, effect_measure)
  }else if(all(names(parameters) %in% c("alpha", "steps"))){
    syntax <-
      "for(j in 1:J){
       eta[j] ~ dgamma(prior_omega_alpha[j], 1)
     }
       for(j in 1:J){
       std_eta[j]  = eta[j] / sum(eta)
       omega[j]    = sum(std_eta[1:j])
     }\n"
  }else if(all(names(parameters) %in% c("alpha1", "alpha2", "steps"))){
    syntax <-
      "for(j1 in 1:J1){
        eta1[j1] ~ dgamma(prior_omega_alpha1[j1], 1)
       }
       for(j2 in 1:J2){
         eta2[j2] ~ dgamma(prior_omega_alpha2[j2], 1)
       }
       for(j1 in 1:J1){
         std_eta1[j1]       = eta1[j1] / sum(eta1)
         omega[J2 - 1 + j1] = sum(std_eta1[1:j1])
       }
         for(j2 in 1:J2){
         std_eta2[j2]  = (eta2[j2] / sum(eta2)) * (1 - std_eta1[1])
       }
       for(j2 in 2:J2){
         omega[j2-1] = sum(std_eta2[j2:J2]) + std_eta1[1]
       }\n"
  }

  return(syntax)
}

# new
.JAGS_survival_mu         <- function(predictors, type){
  mu <- paste0("mu_", type, " = ", "intercept")
  if(!is.null(predictors)){
    mu <- paste0(mu, " + ", paste("beta_", predictors, "*", "x_", type, "_", predictors, collapse = " + ", sep = ""))
  }
  return(paste0(mu, "\n"))
}
.JAGS_survival_likelihood <- function(distribution, type, intercept_only){
  paste0(
    "for(i in 1:", "n_", type, "){\n",
    "   t_", type , "[i]", " ~ ", gsub("-", "_", distribution), "_", type, "(", "mu", "_", type, if(!intercept_only)"[i]",
    if(.has_aux(distribution)){", aux"}, ")\n",
    "}\n")
}
