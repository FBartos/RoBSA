predict.RoBSA        <- function(object, newdata, type = c("survival", "hazard", "density", "mean", "sd"), summarize = TRUE, models = FALSE, conditional = FALSE, samples = 10000){

  if(!(is.matrix(newdata) | is.data.frame(newdata))){
    cnames  <- names(newdata)
    newdata <- data.frame(matrix(newdata, ncol = length(newdata)))
    colnames(newdata) <- cnames
  }
  if(!all(attr(object[["data"]], "predictors") %in% colnames(newdata)))
    stop("All predictors must be provided.")
  if(type %in% c("survival", "hazard", "density")){
    if(!any("time" %in% colnames(newdata)))
      stop("Time variable 'time' must be provided.")
  }


  # deal with predictions for individual models
  if(models){
    out <- lapply(object[["models"]], function(model).predict.RoBSA_model(model, newdata, type, summarize))
    return(out)
  }


  # deal with ensemble level predictions
  out_models <- lapply(object[["models"]], function(model).predict.RoBSA_model(model, newdata, type, FALSE))

  if(conditional){

    predictors    <- attr(object$data, "predictors")
    models        <- object[["models"]]
    marg_liks <- sapply(models, function(x)x$marg_lik$logml)
    mm_predictors <- list()
    for(i in seq_along(predictors)){
      mm_predictors[[predictors[i]]] <- sapply(models, function(m)m$priors[["predictors"]][[predictors[i]]][["type"]] == "alt")
    }
    prior_weights_all        <- sapply(models, function(m)m$prior_odds)
    prior_weights_predictors <- list()
    for(i in seq_along(predictors)){
      prior_weights_predictors[[predictors[i]]] <- ifelse(mm_predictors[[predictors[i]]], prior_weights_all, 0)
    }
    for(i in seq_along(predictors)){
      prior_weights_predictors[[predictors[i]]] <- prior_weights_predictors[[predictors[i]]]/sum(prior_weights_predictors[[predictors[i]]])
    }
    weights_predictors <- list()
    for(i in seq_along(predictors)){
      if(any(mm_predictors[[predictors[i]]]) & all(!is.nan(prior_weights_predictors[[predictors[i]]]))){
        weights_predictors[[predictors[i]]] <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_predictors[[predictors[i]]])
      }
    }

    out <- list()
    for(p in seq_along(predictors)){
      out[[predictors[p]]] <- do.call(cbind, lapply(1:length(weights_predictors[[predictors[p]]]), function(i){
        out_models[[i]][,sample(ncol(out_models[[i]]), round(weights_predictors[[predictors[p]]][i] * samples), TRUE), drop = FALSE]
      }))
    }

    if(summarize){
      for(p in seq_along(predictors)){
        out[[p]] <- data.frame(t(apply(out[[p]], 1, function(x)c(
          "mean" = mean(x),
          "lCI"  = unname(quantile(x, probs = .025)),
          "uCI"  = unname(quantile(x, probs = .975)))
        )))
      }
    }

  }else{

    weights <- object[["RoBSA"]][["posterior_prob"]][["all"]]
    out     <- do.call(cbind, lapply(1:length(weights), function(i){
      out_models[[i]][,sample(ncol(out_models[[i]]), round(weights[i] * samples), TRUE), drop = FALSE]
    }))

    if(summarize){
      out <- data.frame(t(apply(out, 1, function(x)c(
        "mean" = mean(x),
        "lCI"  = unname(quantile(x, probs = .025)),
        "uCI"  = unname(quantile(x, probs = .975)))
      )))
    }

  }


  return(out)
}
.predict.RoBSA_model <- function(model, newdata, type, summarize){

  out <- lapply(1:nrow(newdata), function(i).predict.RoBSA_fun(model, newdata[i,,drop = FALSE], type))
  out <- do.call(rbind, out)

  if(summarize){
    out <- data.frame(t(apply(out, 1, function(x)c(
      "mean" = mean(x),
      "lCI"  = unname(quantile(x, probs = .025)),
      "uCI"  = unname(quantile(x, probs = .975)))
    )))
  }

  return(out)
}
.predict_function    <- function(x, mu, aux, distribution, type){

  # TODO: add other types of censoring
  temp_function <- eval(parse(text = paste0(gsub("-", "_", distribution), "_", type)))

  if(type %in% c("mean", "sd")){
    args <- list(
      eta = mu
    )
  }else{
    args <- list(
      t   = x,
      eta = mu
    )
  }
  if(.has_aux(distribution)){
    args <- c(args, list(aux))
  }

  outcome <- do.call(temp_function, args)

  return(outcome)
}
.predict.RoBSA_fun   <- function(model, data, type){

  priors       <- model[["priors"]]
  samples      <- as.data.frame(suppressWarnings(coda::as.mcmc(model[["fit"]])))
  distribution <- model[["distribution"]]

  if(type %in% c("mean", "sd")){
    x <- NULL
  }else{
    x <- rep(data[["time"]], nrow(samples))
  }

  ### prepare the samples
  # intercept
  if(priors[["intercept"]]$distribution != "point"){
    if(priors[["intercept"]]$distribution == "invgamma"){
      inv_intercept <- samples[,"inv_intercept"]
      intercept     <- 1/inv_intercept
    }else{
      intercept     <- samples[,"intercept"]
    }
  }else{
    intercept <- rep(priors[["intercept"]]$parameters[["location"]], nrow(samples))
  }

  # auxiliary
  if(.has_aux(distribution)){
    if(priors[["aux"]]$distribution != "point"){
      if(priors[["aux"]]$distribution == "invgamma"){
        inv_aux <- samples[,"inv_aux"]
        aux     <- 1/inv_aux
      }else{
        aux     <- samples[,"aux"]
      }
    }else{
      aux <- rep(priors[["aux"]]$parameters[["location"]], nrow(samples))
    }
  }else{
    aux <- NULL
  }

  # predictors
  beta <- list()
  for(i in seq_along(priors[["predictors"]])){
    temp_beta <- NULL

    if(priors[["predictors"]][[i]]$distribution != "point"){
      if(priors[["predictors"]][[i]]$distribution == "invgamma"){
        temp_beta <- samples[,paste0("inv_beta_", names(priors[["predictors"]])[i])]
        temp_beta <- 1/temp_beta
      }else{
        temp_beta <- samples[,paste0( "beta_", names(priors[["predictors"]])[i])]
      }
    }else{
      temp_beta <- rep(priors[["predictors"]][[i]]$parameters[["location"]], nrow(samples))
    }

    beta[[names(priors[["predictors"]])[i]]] <- temp_beta
  }


  # compute the linear predictor
  mu <- intercept
  for(i in seq_along(priors[["predictors"]])){
    mu <- mu + beta[[names(priors[["predictors"]])[i]]] * data[[names(priors[["predictors"]])[i]]]
  }

  out <- .predict_function(x, mu, aux, distribution, type)

  return(out)
}
