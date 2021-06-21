sequential_RoBSA <- function(formula, data, priors = NULL, test_predictors = NULL,
                             distributions = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"),
                             distributions_odds       = rep(1, length(distributions)),
                             default_prior_beta_null  = get_default_prior_beta_null(),
                             default_prior_beta_alt   = get_default_prior_beta_alt(),
                             default_prior_intercept  = get_default_prior_intercept(),
                             default_prior_aux        = get_default_prior_aux(),
                             prior_prob_distributions = rep(1/length(distributions), length(distributions)),
                             nsim = 100000, parallel = FALSE, seed = NULL){

  object       <- NULL
  object$call  <- match.call()

  object$data    <- .prepare_data(formula, data)
  object$priors  <- .prepare_priors(object$priors, distributions, attr(object$data, "predictors"), test_predictors,
                                    default_prior_beta_null, default_prior_beta_alt,
                                    default_prior_intercept, default_prior_aux,
                                    distributions_odds)
  object$models  <- .prepare_models(object$priors)
  object$control <- list(nsim = nsim, parallel = parallel, seed = seed)

  ### fit the models and compute marginal likelihoods
  if(!object$control$parallel){

    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    for(i in 1:length(object$models)){
      object$models[[i]] <- .pointwise_marglik_RoBSA(object, i)
    }

  }else{

    fitting_order <- .fitting_priority(object$models)

    cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
    parallel::clusterEvalQ(cl, {library("RoBSA")})
    parallel::clusterExport(cl, "object", envir = environment())
    object$models <- parallel::clusterApplyLB(cl, fitting_order, .pointwise_marglik_RoBSA, object = object)[order(fitting_order)]
    parallel::stopCluster(cl)

  }

  object$RoBSA         <- .sequential_inference(object)

  class(object) <- "RoBSA.sequential"
  return(object)
}


.pointwise_marglik_RoBSA      <- function(object, i){

  model    <- object$models[[i]]
  priors   <- model$priors
  fit_data <- .fit_data(priors, model$distribution, object$data, model$predictors)

  model$marg_lik <- .pointwise_marglik_function(fit_data, priors, model$distribution, object$control$nsim, object$control$seed)

  return(model)
}
.pointwise_marglik_function   <- function(data, priors, distribution, nsim, seed){

  if(!is.null(seed)){
    set.seed(seed)
  }

  ### simulate from priors
  # intercept and auxiliary parameters
  intercept <- .sim_prior(priors[["intercept"]], nsim)
  aux       <- if(.has_aux(distribution)) .sim_prior(priors[["aux"]], nsim) else NULL

  # predictors
  beta <- list()
  for(i in seq_along(priors[["predictors"]])){
    beta[[names(priors[["predictors"]])[i]]] <- .sim_prior(priors[["predictors"]][[i]], nsim)
  }

  # compute the linear predictor
  mu <- list()
  for(type in attr(data ,"type")){
    mu[[type]] <- .rep_order(intercept, data[[paste0("n_", type)]])
    for(i in seq_along(priors[["predictors"]])){
      mu[[type]] <- mu[[type]] + .rep_order(beta[[names(priors[["predictors"]])[i]]], data[[paste0("n_", type)]]) * rep(data[[paste0("x_", type, "_", names(priors[["predictors"]])[i])]], nsim)
    }
  }

  # compute pointwise marginal log_likelihood
  log_lik <- list()
  for(type in attr(data ,"type")){
    log_lik[[type]] <- .marglik_survival(rep(data[[paste0("t_", type)]], nsim), mu[[type]], .rep_order(aux, data[[paste0("n_", type)]]), distribution, type)
    log_lik[[type]] <- matrix(log_lik[[type]], nrow = data[[paste0("n_", type)]], ncol = nsim)
  }

  # # observation-wise aggregation
  # for(type in attr(data ,"type")){
  #   log_lik[[type]] <- matrix(log_lik[[type]], nrow = data[[paste0("n_", type)]], ncol = nsim)
  #   log_lik[[type]] <- log(apply(exp(log_lik[[type]]), 1, mean))
  # }

  return(log_lik)
}
.pointwise_marglik_sequential <- function(model, data){

  fit_data <- .fit_data(model$priors, model$distribution, data, model$predictors)

  marg_lik <- NULL
  times    <- NULL
  for(type in attr(data ,"type")){
    marg_lik <- rbind(marg_lik, model$marg_lik[[type]])
    times    <- c(times,    data[[paste0("t_", type)]])
  }

  marg_lik <- marg_lik[order(times),]
  times    <- times[order(times)]

  cum_marg_lik <- apply(marg_lik, 2, cumsum)
  cum_marg_lik <- apply(exp(cum_marg_lik), 1, mean)
  cum_marg_lik <- log(cum_marg_lik)

  attr(cum_marg_lik, "times") <- times
  return(cum_marg_lik)
}
.sequential_inference       <- function(object){

  models    <- object$models
  data      <- object$data
  converged <- object$add_info$converged
  seed      <- object$control$seed

  predictors    <- attr(object$data, "predictors")
  distributions <- as.character(sapply(1:length(object$models), function(i)object$models[[i]]$distribution))

  # order the margliks according to the observed time
  cum_marglik <- list()
  for(i in 1:length(models)){
    cum_marglik[[i]] <- .pointwise_marglik_sequential(models[[i]], data)
  }
  ordered_times   <- attr(cum_marglik[[1]], "times")
  cum_marglik <- do.call(cbind, cum_marglik)

  # determine the type of the models
  mm_predictors <- list()
  for(i in seq_along(predictors)){
    mm_predictors[[predictors[i]]] <- sapply(models, function(m)m$priors[["predictors"]][[predictors[i]]][["type"]] == "alt")
  }
  mm_distributions <- list()
  for(i in seq_along(unique(distributions))){
    mm_distributions[[unique(distributions)[i]]] <- sapply(models, function(m)m$distribution == unique(distributions)[i])
  }

  # extract model weights
  prior_weights_all        <- sapply(models, function(m)m$prior_odds)
  prior_weights_predictors <- list()
  for(i in seq_along(predictors)){
    prior_weights_predictors[[predictors[i]]] <- ifelse(mm_predictors[[predictors[i]]], prior_weights_all, 0)
  }
  prior_weights_distributions <- list()
  for(i in seq_along(unique(distributions))){
    prior_weights_distributions[[unique(distributions)[i]]] <- ifelse(mm_distributions[[unique(distributions)[i]]], prior_weights_all, 0)
  }

  # standardize model weights
  prior_weights_all   <- prior_weights_all   / sum(prior_weights_all)
  for(i in seq_along(predictors)){
    prior_weights_predictors[[predictors[i]]] <- prior_weights_predictors[[predictors[i]]]/sum(prior_weights_predictors[[predictors[i]]])
  }
  for(i in seq_along(unique(distributions))){
    prior_weights_distributions[[unique(distributions)[i]]] <- prior_weights_distributions[[unique(distributions)[i]]]/sum(prior_weights_distributions[[unique(distributions)[i]]])
  }

  ### compute model weights
  # overall
  weights_all        <- t(apply(cum_marglik, 1, function(x)
    bridgesampling::post_prob(x, prior_prob = prior_weights_all)
  ))
  weights_predictors <- list()
  for(i in seq_along(predictors)){
    if(any(mm_predictors[[predictors[i]]]) & all(!is.nan(prior_weights_predictors[[predictors[i]]]))){
      weights_predictors[[predictors[i]]] <- t(apply(cum_marglik, 1, function(x)
        bridgesampling::post_prob(x, prior_prob = prior_weights_predictors[[predictors[i]]])
      ))
    }
  }
  weights_distributions <- list()
  for(i in seq_along(unique(distributions))){
    if(any(mm_distributions[[unique(distributions)[i]]]) & all(!is.nan(prior_weights_distributions[[unique(distributions)[i]]]))){
      weights_distributions[[unique(distributions)[i]]] <- t(apply(cum_marglik, 1, function(x)
        bridgesampling::post_prob(x, prior_prob = prior_weights_distributions[[unique(distributions)[i]]])
      ))
    }
  }

  ### compute inclusion BFs
  BF_predictors    <- list()
  for(i in seq_along(predictors)){
    BF_predictors[[predictors[i]]] <- apply(weights_all, 1, function(x)
      .inclusion_BF(prior_weights_all, x, mm_predictors[[predictors[i]]])
    )
  }
  BF_distributions <- list()
  for(i in seq_along(unique(distributions))){
    BF_distributions[[unique(distributions)[i]]] <- apply(weights_all, 1, function(x)
      .inclusion_BF(prior_weights_all, x, mm_distributions[[unique(distributions)[i]]])
    )
  }


  ### edit names
  names(weights_all)   <- names(models)

  ### compute prior and posterior probabilities for model types
  prior_prob_all            <- prior_weights_all
  prior_prob_distributions  <- sapply(unique(distributions), function(distribution)sum(prior_weights_all[mm_distributions[[distribution]]]))
  prior_prob_predictors     <- sapply(predictors,            function(predictor)   sum(prior_weights_all[mm_predictors[[predictor]]]))

  posterior_prob_all            <- weights_all
  posterior_prob_distributions  <- sapply(unique(distributions), function(distribution){
    apply(weights_all, 1, function(x) sum(x[mm_distributions[[distribution]]]))
  })
  posterior_prob_predictors     <- sapply(unique(predictors), function(predictor){
    apply(weights_all, 1, function(x) sum(x[mm_predictors[[predictor]]]))
  })

  output <- list(
    times          = ordered_times,
    BF             = list(
      distributions = BF_distributions,
      predictors    = BF_predictors
    ),
    prior_prob     = list(
      all               = prior_prob_all,
      distributions     = prior_prob_distributions,
      predictors        = prior_prob_predictors
    ),
    posterior_prob = list(
      all               = posterior_prob_all,
      distributions     = posterior_prob_distributions,
      predictors        = posterior_prob_predictors
    )
  )

  return(output)
}


.rep_order   <- function(x, n){
  unlist(lapply(x, function(xi)rep(xi, n)))
}
.sim_prior   <- function(prior, nsim){

  x  <- NULL

  # generate the value
  if(prior$distribution == "normal"){
    while(length(x) < nsim){
      temp_x <- stats::rnorm(nsim, mean = prior$parameters$mean, sd = prior$parameters$sd)
      x <- c(x, temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper])
    }
  }else if(prior$distribution == "lognormal"){
    while(length(x) < nsim){
      temp_x <- stats::rlnorm(nsim, mean = prior$parameters$meanlog, sd = prior$parameters$sdlog)
      x <- c(x, temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper])
    }
  }else if(prior$distribution == "t"){
    while(length(x) < nsim){
      temp_x <- extraDistr::rlst(nsim, df = prior$parameters$df, mu = prior$parameters$location, sigma = prior$parameters$scale)
      x <- c(x, temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper])
    }
  }else if(prior$distribution == "gamma"){
    while(length(x) < nsim){
      temp_x <- stats::rgamma(nsim, shape = prior$parameters$shape, rate = prior$parameters$rate)
      x <- c(x, temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper])
    }
  }else if(prior$distribution == "invgamma"){
    while(length(x) < nsim){
      temp_x <- extraDistr::rinvgamma(nsim, alpha = prior$parameters$shape, beta = prior$parameters$scale)
      x <- c(x, temp_x[temp_x >= prior$truncation$lower & temp_x <= prior$truncation$upper])
    }
  }else if(prior$distribution == "uniform"){
    x <- stats::runif(nsim, min = prior$parameters$a, max = prior$parameters$b)
  }else if(prior$distribution == "point"){
    x <- rep(prior$parameters$location, nsim)
  }

  x <- x[1:nsim]
  return(x)
}
