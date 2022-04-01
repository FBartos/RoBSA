.prepare_data       <- function(formula, data){

  return(list(
    survival   = .prepare_survival(formula, data),
    predictors = .prepare_predictors(formula, data)
  ))
}
.prepare_survival   <- function(formula, data){

  # obtain the survival object
  lhs <- formula[[2]]
  if(lhs[[1]] != "Surv")
    stop("The left hand side of the formula must be a Surv() object.")

  # get the data matrix
  surv_outcome <- eval(lhs, envir = data)

  if(anyNA(surv_outcome))
    stop("NA values are not allowed in the Surv() object.")

  # prepare the data holders
  n_event <- n_cens_r <- n_cens_l <- n_cens_i  <- 0
  t_event <- t_cens_r <- t_cens_l <- t_cens_il <- t_cens_ir <- NULL
  status  <- NULL

  if(attr(surv_outcome, "type") == "right"){

    n_event  <- sum(surv_outcome[,2] == 1)
    n_cens_r <- sum(surv_outcome[,2] == 0)

    t_event  <- surv_outcome[surv_outcome[,2] == 1,1]
    t_cens_r <- surv_outcome[surv_outcome[,2] == 0,1]

  }else if(attr(surv_outcome, "type") == "left"){

    n_event  <- sum(surv_outcome[,2] == 1)
    n_cens_l <- sum(surv_outcome[,2] == 0)

    t_event  <- surv_outcome[surv_outcome[,2] == 1,1]
    t_cens_l <- surv_outcome[surv_outcome[,2] == 0,1]

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

  return(data_survival)
}
.prepare_predictors <- function(formula, data){

  # remove the response
  formula[[2]] <- NULL

  model_frame <- stats::model.frame(formula, data = data)
  model_frame <- as.list(model_frame)
  model_frame

  # add attributes
  data_predictors <- model_frame[1:length(model_frame)]
  attr(data_predictors, "variables") <- attr(attr(model_frame, "terms"), "term.labels")[attr(attr(model_frame, "terms"), "order") == 1]
  attr(data_predictors, "terms")     <- attr(attr(model_frame, "terms"), "term.labels")

  # check for reserved words
  if(any(attr(data_predictors, "terms") %in% c("aux", "intercept", "terms")))
    stop(paste0("The following variable names are internally reserved keywords and cannot be used: ",
                paste0(" '", attr(data_predictors, "terms")[attr(data_predictors, "terms") %in% c("aux", "intercept", "terms")], "' ", collapse = ", ")))

  return(data_predictors)
}
