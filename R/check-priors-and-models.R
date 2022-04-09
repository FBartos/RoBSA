### functions for creating model objects
.check_and_list_priors <- function(priors, distributions, data, test_predictors,
                                   default_prior_beta_null, default_prior_beta_alt,
                                   default_prior_factor_null, default_prior_factor_alt,
                                   default_prior_intercept, default_prior_aux){

  predictors      <- attr(data[["predictors"]],"terms")
  predictors_type <- attr(data[["predictors"]],"terms_type")

  # check the input
  if(!BayesTools::is.prior(default_prior_beta_null))
    stop("The default prior for predictors (null) is not a valid prior distribution.")
  if(!BayesTools::is.prior(default_prior_beta_alt))
    stop("The default prior for predictors (alt) is not a valid prior distribution.")
  if(!all(sapply(default_prior_factor_null, function(p) BayesTools::is.prior.factor(p) | BayesTools::is.prior.point(p))))
    stop("The default prior for factors (null) is not a valid prior distribution.")
  if(!all(sapply(default_prior_factor_alt, function(p) BayesTools::is.prior.factor(p) | BayesTools::is.prior.point(p))))
    stop("The default prior for factors (alt) is not a valid prior distribution.")
  if(!all(sapply(default_prior_intercept[distributions], BayesTools::is.prior)))
    stop("The default prior for intercepts are not a valid prior distribution.")
  if(!all(sapply(default_prior_aux[distributions[.has_aux(distributions)]], BayesTools::is.prior)))
    stop("The default prior for auxilary parameters are not a valid prior distribution.")

  # check for reserved words
  if(any(names(priors) %in% .reserved_words()))
    stop(paste0("The following prior names are internally reserved keywords and cannot be used: ",
                paste0(" '", names(priors)[names(priors) %in% .reserved_words()], "' ", collapse = ", ")))

  # completely the prior distribution specification
  if(is.null(priors) & is.null(test_predictors)){
    # complete default - tests all predictors with default priors
    test_predictors <- predictors

  }else if(!is.null(priors) & is.null(test_predictors)){
    # find whether user specified some parameter priors, if not - tests all predictors with default priors
    test_predictors <- unlist(sapply(predictors, function(p){
      if(!is.null(priors[[p]])){
        if(BayesTools::is.prior(priors[[p]])){
          return("one-prior-is-specified")
        }else{
          if(length(priors[[p]]) == 2 & all(names(priors[[p]]) %in% c("null", "alt"))){
            if(all(sapply(priors[[p]], function(pp) BayesTools::is.prior))){
              return(p)
            }else{
              stop(paste0("The prior distribution for '",p,"' is specified incorrectly."))
            }
          }else{
            stop(paste0("The prior distribution for '",p,"' is specified incorrectly."))
          }
        }
      }else{
        return("no-priors-are-specified")
      }
    }))
    if(all(test_predictors == "no-priors-are-specified")){
      test_predictors <- predictors
    }else{
      test_predictors <- names(predictors[!predictors %in% c("one-prior-is-specified", "no-priors-are-specified")])
    }
  }

  # some additional checks for the remaining specifiable parameters
  if(!is.null(priors)){
    if(!is.null(priors[["aux"]])){
      sapply(distributions, function(d){
        if(!is.null(priors[["aux"]][[d]])){
          if(!BayesTools::is.prior(priors[["aux"]][[d]])){
            stop(paste0("The prior distribution for the auxilary parameter of '",d,"' distribution is specified incorrectly."))
          }
        }
      })
    }
    if(!is.null(priors[["intercept"]])){
      sapply(distributions, function(d){
        if(!is.null(priors[["intercept"]][[d]])){
          if(!BayesTools::is.prior(priors[["intercept"]][[d]])){
            stop(paste0("The prior distribution for the intercept parameter of '",d,"' distribution is specified incorrectly."))
          }
        }
      })
    }
  }


  if(is.null(priors)){

    priors <- list()

    to_test <- predictors[predictors %in% test_predictors]
    no_test <- predictors[!predictors %in% test_predictors]

    for(i in seq_along(to_test)){
      priors[[to_test[i]]] <- list(
        null = if(predictors_type[to_test[i]] == "factor") default_prior_factor_null[["treatment"]] else default_prior_beta_null,
        alt  = if(predictors_type[to_test[i]] == "factor") default_prior_factor_alt[["treatment"]]  else default_prior_beta_alt
      )
    }
    for(i in seq_along(no_test)){
      priors[[no_test[i]]] <- list(
        alt  = if(predictors_type[no_test[i]] == "factor") default_prior_factor_alt[["treatment"]] else default_prior_beta_alt
      )
    }

    priors[["intercept"]] <- default_prior_intercept[distributions]
    priors[["aux"]]       <- default_prior_aux[distributions]

  }else{

    if(any(!names(priors) %in% c(predictors, "intercept"  ,"aux")))
      stop(paste0("The following priors do not corresponds to any predictor or additional parameter: '", paste(names(priors)[!names(priors) %in% c(predictors, "intercept"  ,"aux")], collapse = "', '", sep = ""), "'"))


    to_test <- predictors[predictors %in% test_predictors]
    no_test <- predictors[!predictors %in% test_predictors]

    for(i in seq_along(to_test)){
      if(is.null(priors[[to_test[i]]])){
        priors[[to_test[i]]] <- list(
          null = if(predictors_type[to_test[i]] == "factor") default_prior_factor_null[["treatment"]] else default_prior_beta_null,
          alt  = if(predictors_type[to_test[i]] == "factor") default_prior_factor_alt[["treatment"]]  else default_prior_beta_alt
        )
      }else{
        if(BayesTools::is.prior(priors[[to_test[i]]])){
          priors[[to_test[i]]] <- list(
            null  = if(predictors_type[to_test[i]] == "factor") default_prior_factor_null[[if(BayesTools::is.prior.dummy(priors[[to_test[i]]])) "treatment" else "orthonormal"]] else default_prior_beta_null,
            alt   = priors[[to_test[i]]]
          )
        }else{
          # should be stopped before
          stop(paste0("The predictor '", to_test[i], "' is supposed to be used for testing and the prior distributions are not specified properly"))
        }
      }
    }
    for(i in seq_along(no_test)){
      if(is.null(priors[[no_test[i]]])){
        priors[[no_test[i]]] <- list(
          alt  = default_prior_beta_alt
        )
      }else{
        if(BayesTools::is.prior(priors[[no_test[i]]])){
          priors[[no_test[i]]] <- list(
            alt  = priors[[no_test[i]]]
          )
        }else{
          # should be stopped before
          stop(paste0("The predictor '", no_test[i], "' is supposed to be used for testing and the prior distributions are not specified properly"))
        }
      }
    }
    if(is.null(priors[["intercept"]])){
      priors[["intercept"]] <- default_prior_intercept[distributions]
    }else{
      for(d in distributions){
        if(is.null(priors[["intercept"]][[d]])){
          priors[["intercept"]][[d]] <- default_prior_intercept[d]
        }
      }
    }
    if(is.null(priors[["aux"]])){
      priors[["aux"]] <- default_prior_aux[distributions]
    }else{
      for(d in distributions){
        if(is.null(priors[["aux"]][[d]])){
          priors[["aux"]][[d]] <- default_prior_aux[d]
        }
      }
    }
  }

  priors$terms       <- priors[predictors]
  priors[predictors] <- NULL

  attr(priors, "distributions")  <- distributions
  attr(priors, "terms")          <- predictors
  attr(priors, "terms_test")     <- test_predictors

  return(priors)
}
.prepare_models <- function(priors, distributions, distributions_weights){

  BayesTools::check_char(distributions, "distributions", allow_values = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"), check_length = FALSE)
  BayesTools::check_real(distributions_weights, "distributions_weights", lower = 0, check_length = length(distributions))

  ### create grid of the models
  grid <- list(
    distribution = distributions
  )

  no_test <- attr(priors, "terms")[!attr(priors, "terms") %in% attr(priors, "terms_test")]
  to_test <- attr(priors, "terms")[ attr(priors, "terms") %in% attr(priors, "terms_test")]

  for(i in seq_along(no_test)){
    grid[[no_test[i]]] <- "alt"
  }
  for(i in seq_along(to_test)){
    grid[[to_test[i]]] <- c("null", "alt")
  }

  grid       <- do.call(expand.grid, grid)
  grid$order <- 1:nrow(grid)
  grid       <- merge(grid, data.frame(cbind(distribution = distributions, prior_weights = distributions_weights)), by = "distribution", all.x = TRUE)
  grid       <- grid[order(grid$order),]
  grid$order <- NULL

  if(nrow(grid) > 50)
    warning("More than 50 models are about to be estimated based on the prior specification.", immediate. = TRUE)

  ### create empty models objects for fitting
  models <- lapply(1:nrow(grid), function(i).create_model(grid[i,], priors))

  return(models)
}
.create_model   <- function(grid_row, priors){

  distribution  <- as.character(grid_row[,"distribution"])
  prior_weights <- as.numeric(grid_row[,"prior_weights"])
  terms         <- attr(priors, "terms")

  ### priors
  model_priors <- list()
  model_priors[["intercept"]]  <- priors[["intercept"]][[distribution]]
  if(.has_aux(distribution)){
    model_priors[["aux"]]      <- priors[["aux"]][[distribution]]
  }
  model_priors[["terms"]] <- list()
  terms_test              <- NULL
  for(i in seq_along(terms)){
    model_priors[["terms"]][[terms[i]]]  <- priors[["terms"]][[terms[i]]][[grid_row[,terms[i]]]]
    prior_weights                        <- prior_weights * priors[["terms"]][[terms[i]]][[grid_row[,terms[i]]]]$prior_weights
    if(grid_row[,terms[i]] == "alt"){
      terms_test <- c(terms_test, terms[i])
    }
  }

  model <- list(
    priors       = model_priors,
    distribution = distribution,
    terms        = terms,
    terms_test   = terms_test,
    prior_weights     = prior_weights,
    prior_weights_set = prior_weights
  )
  class(model) <- "RoBSA.model"

  return(model)
}
