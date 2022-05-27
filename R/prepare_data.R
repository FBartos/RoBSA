.prepare_data <- function(formula, data, rescale_data){

  if(!is.language(formula))
    stop("The 'formula' is not specidied as a formula.")
  if(!is.null(data) && !is.data.frame(data))
    stop("'data' must be an object of type data.frame.")
  BayesTools::check_bool(rescale_data, "rescale_data")

  ### obtain the survival object
  lhs <- formula[[2]]
  if(lhs[[1]] != "Surv")
    stop("The left hand side of the formula must be a Surv() object.")

  # get the data matrix
  surv_outcome <- eval(lhs, envir = data)

  if(anyNA(surv_outcome))
    stop("NA values are not allowed in the Surv() object.")

  ### obtain the predictors part
  rhs <- formula[c(1,3)]
  model_frame <- stats::model.frame(rhs, data = data)

  # change characters into factors
  for(i in seq_along(attr(attr(model_frame, "terms"), "dataClasses"))){
    if(attr(attr(model_frame, "terms"), "dataClasses")[[i]] == "character"){
      model_frame[,names(attr(attr(model_frame, "terms"), "dataClasses"))[i]] <-
        as.factor(model_frame[,names(attr(attr(model_frame, "terms"), "dataClasses"))[i]])
      attr(attr(model_frame, "terms"), "dataClasses")[[i]] <- "factor"
    }
  }


  # prepare the data holders
  n_event <- n_cens_r <- n_cens_l <- n_cens_i  <- 0
  t_event <- t_cens_r <- t_cens_l <- t_cens_il <- t_cens_ir <- NULL
  status  <- NULL

  ### sort and order the data
  # always in the following order event, cens_r, cens_l, cens_i
  # (so a single formula with common parameters can be used for model building)
  if(attr(surv_outcome, "type") == "right"){

    n_event  <- sum(surv_outcome[,2] == 1)
    n_cens_r <- sum(surv_outcome[,2] == 0)

    t_event  <- surv_outcome[surv_outcome[,2] == 1,1]
    t_cens_r <- surv_outcome[surv_outcome[,2] == 0,1]

    model_frame <- model_frame[c(
      c(1:nrow(model_frame))[surv_outcome[,2] == 1],
      c(1:nrow(model_frame))[surv_outcome[,2] == 0]
    ),,drop=FALSE]

  }else if(attr(surv_outcome, "type") == "left"){

    n_event  <- sum(surv_outcome[,2] == 1)
    n_cens_l <- sum(surv_outcome[,2] == 0)

    t_event  <- surv_outcome[surv_outcome[,2] == 1,1]
    t_cens_l <- surv_outcome[surv_outcome[,2] == 0,1]

    model_frame <- model_frame[c(
      c(1:nrow(model_frame))[surv_outcome[,2] == 1],
      c(1:nrow(model_frame))[surv_outcome[,2] == 0]
    ),,drop=FALSE]

  }else if(attr(surv_outcome, "type") == "interval"){

    n_event  <- sum(surv_outcome[,3] == 1)
    n_cens_r <- sum(surv_outcome[,3] == 0)
    n_cens_l <- sum(surv_outcome[,3] == 2)
    n_cens_i <- sum(surv_outcome[,3] == 3)

    t_event   <- surv_outcome[surv_outcome[,3] == 1,1]
    t_cens_r  <- surv_outcome[surv_outcome[,3] == 0,1]
    t_cens_l  <- surv_outcome[surv_outcome[,3] == 2,1]
    t_cens_il <- surv_outcome[surv_outcome[,3] == 3,1]
    t_cens_ir <- surv_outcome[surv_outcome[,3] == 3,2]

    model_frame <- model_frame[c(
      c(1:nrow(model_frame))[surv_outcome[,2] == 1],
      c(1:nrow(model_frame))[surv_outcome[,2] == 0],
      c(1:nrow(model_frame))[surv_outcome[,2] == 2],
      c(1:nrow(model_frame))[surv_outcome[,2] == 3]
    ),,drop=FALSE]

  }else{
    stop(paste0("The ", attr(surv_outcome, "type"), " survival type is not implemented."))
  }


  data_survival <- list(
    n_event   = n_event,
    n_cens_r  = n_cens_r,
    n_cens_l  = n_cens_l,
    n_cens_i  = n_cens_i,
    t_event   = t_event,
    t_cens_r  = t_cens_r,
    t_cens_l  = t_cens_l,
    t_cens_il = t_cens_il,
    t_cens_ir = t_cens_ir
  )

  class(data_survival) <- "RoBSA.survdata"
  attr(data_survival, "type") <- c(
    if(n_event > 0)  "event",
    if(n_cens_r > 0) "cens_r",
    if(n_cens_l > 0) "cens_l",
    if(n_cens_i > 0) "cens_i")


  model_frame     <- as.list(model_frame)
  if(length(model_frame) == 0){
    data_predictors <- list()
  }else{
    data_predictors <- model_frame[1:length(model_frame)]
  }
  attr(data_predictors, "variables")  <- attr(attr(model_frame, "terms"), "term.labels")[attr(attr(model_frame, "terms"), "order") == 1]
  attr(data_predictors, "terms")      <- attr(attr(model_frame, "terms"), "term.labels")
  attr(data_predictors, "terms_type") <- attr(attr(model_frame, "terms"), "dataClasses")


  # add additional information about the predictors
  data_predictors_info <- list()
  to_warn              <- NULL
  for(i in seq_along(data_predictors)){
    if(attr(data_predictors, "terms_type")[i] == "numeric"){

      data_predictors_info[[names(data_predictors)[i]]] <- list(
        type = "continuous",
        mean = mean(data_predictors[[names(data_predictors)[i]]]),
        sd   = stats::sd(data_predictors[[names(data_predictors)[i]]])
      )

      if(rescale_data){
        data_predictors[[names(data_predictors)[i]]] <- .pred_scale(data_predictors[[names(data_predictors)[i]]], data_predictors_info[[names(data_predictors)[i]]])
      }else if(RoBSA.get_option("check_scaling") && (abs(mean(data_predictors[[names(data_predictors)[i]]])) > 0.01 | abs(1 - stats::sd(data_predictors[[names(data_predictors)[i]]])) > 0.01)){
        to_warn <- c(to_warn, names(data_predictors)[i])
      }

    }else if(attr(data_predictors, "terms_type")[i] == "factor"){

      data_predictors_info[[names(data_predictors)[i]]] <- list(
        type    = "factor",
        default = levels(data_predictors[[names(data_predictors)[i]]])[1],
        levels  = levels(data_predictors[[names(data_predictors)[i]]])
      )

    }
  }
  attr(data_predictors, "variables_info") <- data_predictors_info

  # create the output object
  output <- list(
    survival   = data_survival,
    predictors = data_predictors
  )

  # throw warnings and errors
  if(length(to_warn) > 0){
    scaling_warning <- paste0("The continuous predictors ", paste0("'", to_warn, "'", collapse = ", "), " are not scaled. Note that extra care need to be taken when specifying prior distributions for unscaled predictors.")
    warning(scaling_warning, immediate. = TRUE, call. = FALSE)
    warning("You can suppress this and following warnings via 'RoBSA.options(check_scaling = FALSE)'. To automatically rescale predictors set 'rescale_data = TRUE'.", immediate. = TRUE, call. = FALSE)
    attr(output, "warnings") <- scaling_warning
  }

  # post processing checks
  if(any(!attr(data_survival, "type") %in% c("event", "cens_r")))
    stop("Only right censored observations are supported.")

  # check for reserved words
  if(any(attr(data_predictors, "terms") %in% .reserved_words()))
    stop(paste0("The following variable names are internally reserved keywords and cannot be used: ",
                paste0(" '", attr(data_predictors, "terms")[attr(data_predictors, "terms") %in% .reserved_words()], "' ", collapse = ", ")))

  return(output)
}


.pred_scale    <- function(x, predictor_info){
  (x - predictor_info[["mean"]]) / predictor_info[["sd"]]
}
.pred_unscale  <- function(x, predictor_info){
  x * predictor_info[["sd"]] + predictor_info[["mean"]]
}
