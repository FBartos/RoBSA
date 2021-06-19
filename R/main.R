RoBSA <- function(formula, data, priors = NULL, test_predictors = NULL,
                  distributions = c("exp-aft", "weibull-aft", "lnorm-aft", "llogis-aft", "gamma-aft"),
                  distributions_odds       = rep(1, length(distributions)),
                  default_prior_beta_null  = get_default_prior_beta_null(),
                  default_prior_beta_alt   = get_default_prior_beta_alt(),
                  default_prior_intercept  = get_default_prior_intercept(),
                  default_prior_aux        = get_default_prior_aux(),
                  prior_prob_distributions = rep(1/length(distributions), length(distributions)),
                  chains = 2, iter = 5000, burnin = 1000, thin = 1, parallel = FALSE,
                  control = NULL, save = "all", seed = NULL){

  object       <- NULL
  object$call  <- match.call()

  object$data    <- .prepare_data(formula, data)
  object$priors  <- .prepare_priors(priors, distributions, attr(object$data, "predictors"), test_predictors,
                                   default_prior_beta_null, default_prior_beta_alt,
                                   default_prior_intercept, default_prior_aux,
                                   distributions_odds)
  object$models  <- .prepare_models(object$priors)
  object$control <- .set_control(control, chains, iter, burnin, thin, seed, parallel)

  ### fit the models and compute marginal likelihoods
  if(object$control$cores < 2*object$control$chains){

    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    for(i in 1:length(object$models)){
      object$models[[i]] <- .fit_RoBSA_wrap(object, i)
    }

  }else{

    fitting_order <- .fitting_priority(object$models)

    cl <- parallel::makePSOCKcluster(floor(object$control$cores / object$control$chains))
    parallel::clusterEvalQ(cl, {library("RoBSA")})
    parallel::clusterExport(cl, "object", envir = environment())
    object$models <- parallel::clusterApplyLB(cl, fitting_order, .fit_RoBSA_wrap, object = object)[order(fitting_order)]
    parallel::stopCluster(cl)

  }


    # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)


  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){

    # TODO: any idea how to proceed with this?
    # balance probability of non-converged models TODO
    # if(object$control$balance_prob & any(!object$add_info$converged))object <- .balance_prob(object, object$add_info$converged)


    ### compute the model-space results
    object$RoBSA         <- .model_inference(object)
    object$coefficients  <- .compute_coeficients(object$RoBSA)
  }


  ### add warnings
  object$add_info$warnings <- c(object$add_info$warnings, .model_refit_warnings(sapply(1:length(object$models), function(i)object$models[[i]]$metadata, simplify = FALSE)))
  object$add_info$warnings <- c(object$add_info$warnings, .model_convergence_warnings(object))


  ### remove model posteriors if asked to
  if(save == "min"){
    for(i in 1:length(object$models)){
      if(length(object$models[[1]]$fit) != 1){
        object$models[[i]]$fit$mcmc <- NULL
      }
    }
  }


  ### print warnings
  if(!is.null(object$add_info$warnings)){
    for(w in object$add_info$warnings)
      warning(w)
  }
  if(sum(!object$add_info$converged) > 0)
    warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))


  class(object) <- "RoBSA"
  return(object)
}


update.RoBSA <- function(object, refit_failed = TRUE,
                         priors = NULL, prior_odds = NULL,
                         control = NULL, chains = NULL, iter = NULL, burnin = NULL, thin = NULL, parallel = NULL, seed = NULL, ...){


  # if(object$add_info$save == "min")
  #   stop("Models cannot be updated because individual model posteriors were not save during the fitting process. Set 'save' parameter to 'all' in while fitting the model (see ?RoBSA for more details).")



  ### choose proper action based on the supplied input
  if(!is.null(priors)){

    what_to_do <- "fit_new_model"

    if(!is.null(prior_odds))object$models[[length(object$models)]]$prior_odds     <- prior_odds
    if(!is.null(prior_odds))object$models[[length(object$models)]]$prior_odds_set <- prior_odds


  }else if(!is.null(prior_odds)){

    what_to_do <- "update_prior_odds"
    if(length(prior_odds) != length(object$models))
      stop("The number of newly specified prior odds does not match the number of models. See '?update.RoSMA' for more details.")
    for(i in 1:length(object$models)){
      object$models[[i]]$prior_odds     <- prior_odds[i]
      object$models[[i]]$prior_odds_set <- prior_odds[i]
    }

  }else if(refit_failed & any(!object$add_info$converged)){

    what_to_do <- "refit_failed_models"

  }else{

    what_to_do <- "update_settings"

  }


  ### update control settings if any change is specified
  object$control  <- .update_control(object$control, control, chains, iter, burnin, thin, seed, parallel)


  ### do the stuff
  if(what_to_do == "fit_new_model"){

    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    object$models[[length(object$models)]] <- .fit_RoBSA_wrap(object, length(object$models))

  }else if(what_to_do == "refit_failed_models"){

    converged_models <- .get_converged_models(object)
    if(!is.null(object$control$progress_start))eval(parse(text = object$control$progress_start))
    for(i in c(1:length(object$models))[!converged_models]){
      object$models[[i]] <- .fit_RoBSA_wrap(object, i)
    }

  }else if(what_to_do == "transform_estimates"){

    for(i in c(1:length(object$models))){
      object$models[[i]] <- .transform_posteriors(object$models[[i]], object$add_info$output_scale, .transformation_var(output_scale))
    }
    object <- .transform_posterior(object, object$add_info$output_scale, .transformation_var(output_scale))
    object$add_info$output_scale <- .transformation_var(output_scale)

    return(object)

  }


  # deal with non-converged the converged models
  object$add_info$converged <- .get_converged_models(object)

  # create ensemble only if at least one model converges
  if(any(object$add_info$converged)){

    ### compute the model-space results
    object$RoBSA         <- .model_inference(object)
    object$coefficients  <- .compute_coeficients(object$RoBSA)
  }


  ### add warnings
  object$add_info$warnings <- c(object$add_info$warnings, .model_refit_warnings(sapply(1:length(object$models), function(i)object$models[[i]]$metadata, simplify = FALSE)))
  object$add_info$warnings <- c(object$add_info$warnings, .model_convergence_warnings(object))


  ### remove model posteriors if asked to
  # if(save == "min"){
  #   for(i in 1:length(object$models)){
  #     if(length(object$models[[1]]$fit) != 1){
  #       object$models[[i]]$fit$mcmc <- NULL
  #     }
  #   }
  # }


  ### print warnings
  if(!is.null(object$add_info$warnings)){
    for(w in object$add_info$warnings)
      warning(w)
  }
  if(sum(!object$add_info$converged) > 0)
    warning(paste0(sum(!object$add_info$converged), ifelse(sum(!object$add_info$converged) == 1, " model", " models"), " failed to converge."))

  return(object)
}


.prepare_data   <- function(formula, data){
  ### adapted from rstanarm stan_surv

  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function.")

  ### process data
  if (missing(data) || !inherits(data, "data.frame"))
    stop("'data' must be a data frame.")


  nlist <- function(...) {
    m <- match.call()
    out <- list(...)
    no_names <- is.null(names(out))
    has_name <- if (no_names) FALSE else nzchar(names(out))
    if (all(has_name))
      return(out)
    nms <- as.character(m)[-1L]
    if (no_names) {
      names(out) <- nms
    } else {
      names(out)[!has_name] <- nms[!has_name]
    }

    return(out)
  }
  parse_formula_and_data <- function(formula, data) {

    validate_formula <- function(formula, needs_response = TRUE) {

      if (!inherits(formula, "formula")) {
        stop2("'formula' must be a formula.")
      }

      if (needs_response) {
        len <- length(formula)
        if (len < 3) {
          stop2("'formula' must contain a response.")
        }
      }
      as.formula(formula)
    }
    formula <- validate_formula(formula, needs_response = TRUE)

    # all variables of entire formula
    allvars <- all.vars(formula)
    allvars_form <- reformulate(allvars)

    # LHS of entire formula
    lhs <- function(x, as_formula = FALSE) {
      len <- length(x)
      if (len == 3L) {
        out <- x[[2L]]
      } else {
        out <- NULL
      }
      out
    }
    reformulate_lhs <- function(x) {
      x <- formula(substitute(LHS ~ 1, list(LHS = x)))
      x
    }
    lhs       <- lhs(formula)         # LHS as expression
    lhs_form  <- reformulate_lhs(lhs) # LHS as formula

    # RHS of entire formula
    rhs <- function(x, as_formula = FALSE) {
      len <- length(x)
      if (len == 3L) {
        out <- x[[3L]]
      } else {
        out <- x[[2L]]
      }
      out
    }
    reformulate_rhs <- function(x) {
      x <- formula(substitute(~ RHS, list(RHS = x)))
      x
    }
    rhs       <- rhs(formula)         # RHS as expression
    rhs_form  <- reformulate_rhs(rhs) # RHS as formula

    # evaluate model data (row subsetting etc)
    make_model_data <- function(formula, data) {
      mf <- model.frame(formula, data, na.action = na.pass)
      include <- apply(mf, 1L, function(row) !any(is.na(row)))
      data[include, , drop = FALSE]
    }
    data <- make_model_data(allvars_form, data)

    # evaluated response variables
    validate_surv <- function(x, ok_types = c("right", "counting", "interval", "interval2")) {
      if (!inherits(x, "Surv"))
        stop2("LHS of 'formula' must be a 'Surv' object.")
      if (!attr(x, "type") %in% ok_types)
        stop2("Surv object type must be one of: ", comma(ok_types))
      x
    }

    surv <- eval(lhs, envir = data) # Surv object
    surv <- validate_surv(surv)
    type <- attr(surv, "type")

    status   <- as.vector(surv[, "status"])


    if (any(is.na(status)))
      stop2("Invalid status indicator in Surv object.")

    if (any(status < 0 || status > 3))
      stop2("Invalid status indicator in Surv object.")


    # just fixed-effect part of formula
    fe_form   <- lme4::nobars(formula)

    # just random-effect part of formula
    split_at_bars <- function(x) {
      terms <- strsplit(deparse(x, 500), "\\s\\|\\s")[[1L]]
      if (!length(terms) == 2L)
        stop2("Could not parse the random effects formula.")
      re_form <- formula(paste("~", terms[[1L]]))
      group_var <- terms[[2L]]
      nlist(re_form, group_var)
    }
    fetch <- function(x, y, z = NULL, zz = NULL, null_to_zero = FALSE,
                      pad_length = NULL, unlist = FALSE) {
      ret <- lapply(x, `[[`, y)
      if (!is.null(z))
        ret <- lapply(ret, `[[`, z)
      if (!is.null(zz))
        ret <- lapply(ret, `[[`, zz)
      if (null_to_zero)
        ret <- lapply(ret, function(i) ifelse(is.null(i), 0L, i))
      if (!is.null(pad_length)) {
        padding <- rep(list(0L), pad_length - length(ret))
        ret <- c(ret, padding)
      }
      if (unlist) unlist(ret) else ret
    }
    bars      <- lme4::findbars(formula)
    re_parts  <- lapply(bars, split_at_bars)
    re_forms  <- fetch(re_parts, "re_form")

    tf_form <- formula

    nlist(formula,
          data,
          allvars,
          allvars_form,
          lhs,
          lhs_form,
          rhs,
          rhs_form,
          formula,
          tf_form,
          fe_form,
          bars,
          re_parts,
          re_forms,
          surv_type = attr(surv, "type"))
  }

  formula <- parse_formula_and_data(formula, data)
  data    <- formula$data; formula[["data"]] <- NULL

  #----------------
  # Construct data
  #----------------

  #----- model frame stuff
  make_model_frame <- function(formula,
                               data,
                               xlevs              = NULL,
                               drop.unused.levels = FALSE,
                               check_constant     = FALSE,
                               na.action          = na.fail) {

    # construct model frame
    Terms <- terms(lme4::subbars(formula))
    mf <- stats::model.frame(Terms,
                             data,
                             xlev = xlevs,
                             drop.unused.levels = drop.unused.levels,
                             na.action = na.action)

    # get predvars for fixed part of formula
    TermsF <- terms(lme4::nobars(formula))
    mfF <- stats::model.frame(TermsF,
                              data,
                              xlev = xlevs,
                              drop.unused.levels = drop.unused.levels,
                              na.action = na.action)
    attr(attr(mf, "terms"), "predvars.fixed") <- attr(attr(mfF, "terms"), "predvars")

    # get predvars for random part of formula
    has_bars <- length(lme4::findbars(formula)) > 0
    if (has_bars) {
      TermsR <- terms(lme4::subbars(justRE(formula, response = TRUE)))
      mfR <- stats::model.frame(TermsR,
                                data,
                                xlev = xlevs,
                                drop.unused.levels = drop.unused.levels,
                                na.action = na.action)
      attr(attr(mf, "terms"), "predvars.random") <- attr(attr(mfR, "terms"), "predvars")
    } else {
      attr(attr(mf, "terms"), "predvars.random") <- NULL
    }

    # check no constant vars
    check_constant_vars <- function(mf) {
      mf1 <- mf
      if (NCOL(mf[, 1]) == 2 || all(mf[, 1] %in% c(0, 1)) || survival::is.Surv(mf[, 1])) {
        mf1 <- mf[, -1, drop=FALSE]
      }

      lu1 <- function(x) !all(x == 1) && length(unique(x)) == 1
      nocheck <- c("(weights)", "(offset)", "(Intercept)")
      sel <- !colnames(mf1) %in% nocheck
      is_constant <- apply(mf1[, sel, drop=FALSE], 2, lu1)
      if (any(is_constant)) {
        stop("Constant variable(s) found: ",
             paste(names(is_constant)[is_constant], collapse = ", "),
             call. = FALSE)
      }
      return(mf)
    }
    if (check_constant)
      mf <- check_constant_vars(mf)

    # add additional predvars attributes

    # check for terms
    mt <- attr(mf, "terms")
    if (is.empty.model(mt))
      stop2("No intercept or predictors specified.")

    nlist(mf, mt)
  }
  drop_intercept <- function(x) {
    sel <- check_for_intercept(x)
    if (length(sel) && is.matrix(x)) {
      x[, -sel, drop = FALSE]
    } else if (length(sel)) {
      x[-sel]
    } else {
      x
    }
  }
  check_for_intercept  <- function(x, logical = FALSE) {
    nms <- if (is.matrix(x)) colnames(x) else names(x)
    sel <- which("(Intercept)" %in% nms)
    if (logical) as.logical(length(sel)) else sel
  }

  mf_stuff <- make_model_frame(formula$tf_form, data, drop.unused.levels = TRUE)

  mf <- mf_stuff$mf # model frame
  mt <- mf_stuff$mt # model terms

  #----- dimensions and response vectors

  # entry and exit times for each row of data
  make_t <- function(model_frame, type = c("beg", "end", "gap", "upp")) {

    type <- match.arg(type)
    resp <- if (survival::is.Surv(model_frame))
      model_frame else model.response(model_frame)
    surv <- attr(resp, "type")
    err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

    t_beg <- switch(surv,
                    "right"     = rep(0, nrow(model_frame)),
                    "interval"  = rep(0, nrow(model_frame)),
                    "interval2" = rep(0, nrow(model_frame)),
                    "counting"  = as.vector(resp[, "start"]),
                    stop(err))

    t_end <- switch(surv,
                    "right"     = as.vector(resp[, "time"]),
                    "interval"  = as.vector(resp[, "time1"]),
                    "interval2" = as.vector(resp[, "time1"]),
                    "counting"  = as.vector(resp[, "stop"]),
                    stop(err))

    t_upp <- switch(surv,
                    "right"     = rep(NaN, nrow(model_frame)),
                    "counting"  = rep(NaN, nrow(model_frame)),
                    "interval"  = as.vector(resp[, "time2"]),
                    "interval2" = as.vector(resp[, "time2"]),
                    stop(err))

    switch(type,
           "beg" = t_beg,
           "end" = t_end,
           "gap" = t_end - t_beg,
           "upp" = t_upp,
           stop("Bug found: cannot handle specified 'type'."))
  }
  t_beg <- make_t(mf, type = "beg") # entry time
  t_end <- make_t(mf, type = "end") # exit  time
  t_upp <- make_t(mf, type = "upp") # upper time for interval censoring

  # ensure no event or censoring times are zero (leads to degenerate
  # estimate for log hazard for most baseline hazards, due to log(0))
  check1 <- any(t_end <= 0, na.rm = TRUE)
  check2 <- any(t_upp <= 0, na.rm = TRUE)
  if (check1 || check2)
    stop2("All event and censoring times must be greater than 0.")

  # event indicator for each row of data
  make_d <- function(model_frame) {

    resp <- if (survival::is.Surv(model_frame))
      model_frame else model.response(model_frame)
    surv <- attr(resp, "type")
    err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

    switch(surv,
           "right"     = as.vector(resp[, "status"]),
           "interval"  = as.vector(resp[, "status"]),
           "interval2" = as.vector(resp[, "status"]),
           "counting"  = as.vector(resp[, "status"]),
           stop(err))
  }
  status <- make_d(mf)


  if (any(is.na(status)))
    stop2("Invalid status indicator in Surv object.")

  if (any(status < 0 || status > 3))
    stop2("Invalid status indicator in Surv object.")

  # delayed entry indicator for each row of data
  delayed  <- as.logical(!t_beg == 0)

  # time variables for stan
  t_event <- t_end[status == 1] # exact event time
  t_lcent <- t_end[status == 2] # left  censoring time
  t_rcent <- t_end[status == 0] # right censoring time
  t_icenl <- t_end[status == 3] # lower limit of interval censoring time
  t_icenu <- t_upp[status == 3] # upper limit of interval censoring time
  t_delay <- t_beg[delayed]     # delayed entry time

  # dimensions
  n_event <- sum(status == 1)
  n_rcent <- sum(status == 0)
  n_lcent <- sum(status == 2)
  n_icent <- sum(status == 3)
  n_delay <- sum(delayed)

  # predictors
  make_x <- function(formula,
                     model_frame,
                     xlevs = NULL,
                     check_constant = TRUE) {

    # uncentred predictor matrix, without intercept
    x <- model.matrix(formula, model_frame, xlev = xlevs)
    x <- drop_intercept(x)

    # column means of predictor matrix
    x_bar <- as.array(colMeans(x))

    # centered predictor matrix
    x_centered <- sweep(x, 2, x_bar, FUN = "-")

    # identify any column of x with < 2 unique values (empty interaction levels)
    sel <- (apply(x, 2L, n_distinct) < 2)
    if (check_constant && any(sel)) {
      cols <- paste(colnames(x)[sel], collapse = ", ")
      stop2("Cannot deal with empty interaction levels found in columns: ", cols)
    }

    nlist(x, x_centered, x_bar, N = NROW(x), K = NCOL(x))
  }
  n_distinct  <-  function(x) {
    length(unique(x))
  }

  rep_rows <- function(x, ...) {
    if (is.null(x) || !nrow(x)) {
      return(x)
    } else if (is.matrix(x) || is.data.frame(x)) {
      x <- x[rep(1:nrow(x), ...), , drop = FALSE]
    } else {
      stop2("'x' must be a matrix or data frame.")
    }
    x
  }
  keep_rows <-  function(x, rows = 1:nrow(x)) {
    x[rows, , drop = FALSE]
  }

  mf_event <- keep_rows(mf, status == 1)
  mf_lcent <- keep_rows(mf, status == 2)
  mf_rcent <- keep_rows(mf, status == 0)
  mf_icent <- keep_rows(mf, status == 3)
  mf_delay <- keep_rows(mf, delayed)

  mf_cpts <- rbind(mf_event,
                   mf_lcent,
                   mf_rcent,
                   mf_icent,
                   mf_icent,
                   mf_delay)

  ff        <- formula$fe_form
  x         <- make_x(ff, mf     )$x
  x_cpts    <- make_x(ff, mf_cpts)$x
  x_centred <- sweep(x_cpts, 2, colMeans(x), FUN = "-")
  K         <- ncol(x_cpts)



  add_x <- function(out, x, status, type){
    for(i in seq_along(colnames(x)))
      out[[paste0("x", "_", type, "_", colnames(x)[i])]] <- x[status, colnames(x)[i]]
    out
  }

  out <- list()
  attr(out, "predictors") <- colnames(x)
  if(any(colnames(x) %in% c("intercept", "aux")))
    stop("'intercept' and 'aux' are reserved names and cannot be used as variable names.")

  # right censoring time
  if(n_rcent > 0) {
    out$n_rcent  <- n_rcent
    out$t_rcent <- t_rcent
    out <- add_x(out, x, status == 0, "rcent")
  }
  # exact event time
  if(n_event > 0) {
    out$n_event  <- n_event
    out$t_event <- t_event
    out <- add_x(out, x, status == 1, "event")
  }
  # left  censoring time
  if(n_lcent > 0) {
    out$n_lcent  <- n_lcent
    out$t_lcent <- t_lcent
    out <- add_x(out, x, status == 2, "lcent")
  }
  # lower and upper limit of interval censoring time
  if(n_icent > 0) {
    out$n_icent  <- n_icent
    out$t_icenl <- t_icenl
    out$t_icenu <- t_icenu
    out <- add_x(out, x, status == 3, "icent")
  }
  # delayed entry time
  if(n_delay > 0) {
    out$n_delay  <- n_delay
    out$t_delay <- t_delay
    out <- add_x(out, x, delayed, "delayed")
  }

  attr(out, "type") <- c(
    if(n_event > 0) "event",
    if(n_rcent > 0) "rcent",
    if(n_lcent > 0) "lcent",
    if(n_icent > 0) "icent",
    if(n_delay > 0) "delay"
  )
  return(out)
}
.prepare_priors <- function(priors, distributions, predictors, test_predictors,
                           default_prior_beta_null, default_prior_beta_alt,
                           default_prior_intercept, default_prior_aux,
                           distributions_odds){


  if(!inherits(default_prior_beta_null, "RoBSA.prior")){
    stop("The default prior for predictors (null) is not a valid prior distribution.")
  }
  if(!inherits(default_prior_beta_alt, "RoBSA.prior")){
    stop("The default prior for predictors (alt) is not a valid prior distribution.")
  }
  if(!all(sapply(default_prior_intercept[distributions], function(p)inherits(p, "RoBSA.prior")))){
    stop("The default prior for intercepts are not a valid prior distribution.")
  }
  if(!all(sapply(default_prior_aux[distributions[.has_aux(distributions)]], function(p)inherits(p, "RoBSA.prior")))){
    stop("The default prior for auxilary parameters are not a valid prior distribution.")
  }


  if(is.null(priors) & is.null(test_predictors)){
    # complete default - tests all predictors with default priors
    test_predictors <- predictors

  }else if(!is.null(priors) & is.null(test_predictors)){
    # find whether user specified some parameter priors, if not - tests all predictors with default priors
    test_predictors <- unlist(sapply(predictors, function(p){
      if(!is.null(priors[[p]])){
        if(inherits(priors[[p]], "RoBSA.prior")){
          return("one-prior-is-specified")
        }else{
          if(length(priors[[p]]) == 2 & all(names(priors[[p]]) %in% c("null", "alt"))){
            if(all(sapply(priors[[p]], function(pp)inherits("RoBSA.prior")))){
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
          if(!inherits(priors[[d]][["aux"]], "RoBSA.prior")){
            stop(paste0("The prior distribution for the auxilary parameter of '",d,"' distribution is specified incorrectly."))
          }
        }
      })
    }
    if(!is.null(priors[["intercept"]])){
      sapply(distributions, function(d){
        if(!is.null(priors[["intercept"]][[d]])){
          if(!inherits(priors[[d]][["intercept"]], "RoBSA.prior")){
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
        null = default_prior_beta_null,
        alt  = default_prior_beta_alt
      )
    }
    for(i in seq_along(no_test)){
      priors[[no_test[i]]] <- list(
        alt  = default_prior_beta_alt
      )
    }

    priors[["intercept"]] <- default_prior_intercept[distributions]
    priors[["aux"]]       <- default_prior_aux[distributions]

  }else{

    if(any(!names(priors) %in% c(predictors, "intercept"  ,"aux")))
      stop(paste0("The following priors do not corresponds to any predictor or additional parameter: '", paste(!names(priors) %in% c(predictors, "intercept"  ,"aux"), collapse = "', '", sep = ""), "'"))

    to_test <- predictors[predictors %in% test_predictors]
    no_test <- predictors[!predictors %in% test_predictors]

    for(i in seq_along(to_test)){
      if(is.null(priors[[to_test[i]]])){
        priors[[to_test[i]]] <- list(
          null = default_prior_beta_null,
          alt  = default_prior_beta_alt
        )
      }else{
        if(inherits(priors[[to_test[i]]], "RoBSA.prior")){
          priors[[to_test[i]]] <- list(
            null = default_prior_beta_null,
            alt  = priors[[to_test[i]]]
          )
        }else{
          # should be stoppe before
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
        if(inherits(priors[[no_test[i]]], "RoBSA.prior")){
          priors[[no_test[i]]] <- list(
            alt  = priors[[no_test[i]]]
          )
        }else{
          # should be stoppe before
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

  if(is.null(names(distributions_odds))){
    names(distributions_odds) <- distributions
  }

  attr(priors, "distributions")      <- distributions
  attr(priors, "distributions_odds") <- distributions_odds
  attr(priors, "predictors")         <- predictors
  attr(priors, "test_predictors")    <- test_predictors

  return(priors)
}
.prepare_models <- function(priors){

  ### create grid of the models
  grid <- list(
    distribution = attr(priors, "distributions")
  )

  no_test <- attr(priors, "predictors")[!attr(priors, "predictors") %in% attr(priors, "test_predictors")]
  to_test <- attr(priors, "predictors")[attr(priors, "predictors") %in% attr(priors, "test_predictors")]

  for(i in seq_along(no_test)){
    grid[[no_test[i]]] <- "alt"
  }
  for(i in seq_along(to_test)){
    grid[[to_test[i]]] <- c("null", "alt")
  }

  grid       <- do.call(expand.grid, grid)
  grid$order <- 1:nrow(grid)
  prior_odds <- data.frame(cbind(distribution = names(attr(priors, "distributions_odds")), prior_odds = attr(priors, "distributions_odds")))
  grid       <- merge(grid, prior_odds, by = "distribution", all.x = TRUE)
  grid       <- grid[order(grid$order),]
  grid$order <- NULL

  ### create empty models objects for fitting
  models <- lapply(1:nrow(grid), function(i).create_model(grid[i,], priors))

  return(models)
}
.create_model   <- function(grid_row, priors){

  distribution <- grid_row[,"distribution"]
  prior_odds   <- as.numeric(grid_row[,"prior_odds"])
  predictors   <- attr(priors, "predictors")

  ### priors
  model_priors <- list()
  model_priors[["intercept"]]  <- priors[["intercept"]][[distribution]]
  if(.has_aux(distribution)){
    model_priors[["aux"]]      <- priors[["aux"]][[distribution]]
  }
  model_priors[["predictors"]] <- list()
  for(i in seq_along(predictors)){
    model_priors[["predictors"]][[predictors[i]]] <- priors[[predictors[i]]][[grid_row[,predictors[i]]]]
    model_priors[["predictors"]][[predictors[i]]][["type"]] <- grid_row[,predictors[i]]
    prior_odds <- prior_odds * priors[[predictors[i]]][[grid_row[,predictors[i]]]]$prior_odds
  }

  model <- list(
    priors       = model_priors,
    distribution = distribution,
    predictors   = predictors,
    prior_odds   = prior_odds
  )

  return(model)
}
.fit_RoBSA_wrap <- function(object, i){

  object$models[[i]] <- .fit_RoBSA(object, i)
  object$models[[i]] <- .marglik_RoBSA(object, i)
  if(!is.null(object$control$progress_tick))eval(parse(text = object$control$progress_tick))

  return(object$models[[i]])
}
.fit_RoBSA      <- function(object, i){

  model      <- object$models[[i]]
  priors     <- model$priors
  control    <- object$control

  refit_info <- NULL
  fit        <- NULL


  # generate additional information for model fitting
  model_syntax      <- .generate_model_syntax(priors, model$distribution, object$data)
  fit_data          <- .fit_data(priors, model$distribution, object$data, attr(model_syntax, "predictors"))
  fit_inits         <- .fit_inits(priors, model$distribution, control$chains, control$seed)
  monitor_variables <- .to_monitor(model$distribution, attr(model_syntax, "predictors"))

  # fit the model
  fit <- .fit_model_RoBSA_wrap(model_syntax, fit_data, fit_inits, monitor_variables, control)

  # error handling
  if(all(class(fit) %in% c("simpleError", "error", "condition"))){

    # problem with installing the RoBSA JAGS module
    if(grepl("Unknown distribution", fit$message)){
      stop("The RoBSA JAGS distributions could not be found. Please, check that the RoBSA package is properly installed.")
    }

    # are there any fixable errors?

  }

  # forward error if it's unfixable
  if(all(class(fit) %in% c("simpleError", "error", "condition"))){
    refit_info <- fit$message
  }

  # add the fit and summary to the main object
  model$fit      <- fit
  model$metadata <- list(
    i          = i,
    refit_info = refit_info)
  if(!is.null(fit) & !any(class(fit) %in% c("simpleError", "error", "condition"))){
    model$fit_summary <- .runjags.summary(fit, model$distribution)
  }

  return(model)
}
.marglik_RoBSA  <- function(object, i){

  model    <- object$models[[i]]
  fit      <- model$fit
  priors   <- model$priors
  control  <- object$control

  # deal with failed model
  if(any(class(fit) %in% c("simpleError", "error"))){
    model$marg_lik <- .marglik_fail()
    return(model)
  }

  # compute marginal likelihood
  fit_data        <- .fit_data(priors, model$distribution, object$data, model$predictors)
  marglik_samples <- .marglik_prepare_data(fit, priors, model$distribution, fit_data)

  if(!is.null(control$seed))set.seed(control$seed)
  marg_lik        <- tryCatch(suppressWarnings(bridgesampling::bridge_sampler(
    samples          = marglik_samples$samples,
    data             = fit_data,
    log_posterior    = .marglik_function,
    distribution     = model$distribution,
    priors           = priors,
    lb               = marglik_samples$lb,
    ub               = marglik_samples$ub,
    maxiter          = control$bridge_max_iter,
    silent           = TRUE)),
    error = function(e)return(e))

  # handle errors
  if(any(class(marg_lik) %in% c("simpleError", "error"))){
    model$metadata$marg_lik <- marg_lik$message
    marg_lik <- .marglik_fail()
  }else if(is.na(marg_lik$logml)){
    model$metadata$marg_lik <- "not enough iterations"
    marg_lik <- .marglik_fail()
  }

  model$marg_lik <- marg_lik

  return(model)
}

### model syntax
.generate_model_syntax <- function(priors, distribution, data){

  # generate model syntax
  model_syntax <- "model{\n"

  ### prior distributions
  # intercept and auxiliary parameters
  model_syntax <- paste0(model_syntax, .JAGS_distribution("intercept", priors[["intercept"]]$distribution, priors[["intercept"]]$truncation))
  if(.has_aux(distribution)){
    model_syntax <- paste0(model_syntax, .JAGS_distribution("aux", priors[["aux"]]$distribution, priors[["aux"]]$truncation))
  }

  # predictors
  predictors <- NULL
  for(i in seq_along(priors$predictors)){
    if(!.is_zero_spike(priors$predictors[[i]])){
      predictors   <- c(predictors, names(priors$predictors)[i])
      model_syntax <- paste0(model_syntax, .JAGS_distribution(paste0("beta_", names(priors[["predictors"]])[i]), priors[["predictors"]][[i]]$distribution, priors[["predictors"]][[i]]$truncation))
    }
  }

  ### model
  for(type in attr(data ,"type")){
    if(!is.null(data[[paste0("n_", type)]])){
      model_syntax <- paste0(model_syntax, .JAGS_survival_mu(predictors, type))
    }
  }

  ### likelihoods for the observed data
  for(type in attr(data ,"type")){
    if(!is.null(data[[paste0("n_", type)]])){
      model_syntax <- paste0(model_syntax, .JAGS_survival_likelihood(distribution, type, is.null(predictors)))
    }
  }

  model_syntax <- paste0(model_syntax, "}")

  attr(model_syntax, "predictors") <- predictors

  return(model_syntax)
}
.fit_data              <- function(priors, distribution, data, predictors){

  # remove predictors that are not used
  not_used <- attr(data, "predictors")[!attr(data, "predictors") %in% predictors]
  for(i in seq_along(not_used)){
    for(type in attr(data, "type")){
      data[[paste0("x_",type,"_",not_used[i])]] <- NULL
    }
  }

  ### add prior parameter values
  # intercept and auxiliary parameters
  for(var in c("intercept", if(.has_aux(distribution)) "aux")){
    for(par in names(priors[[var]]$parameters)){
      data[[paste0("prior_",var,"_",par)]] <- priors[[var]]$parameters[[par]]
    }
  }

  # predictors
  for(i in seq_along(priors[["predictors"]])){
    if(!.is_zero_spike(priors[["predictors"]][[i]])){
      for(par in names(priors[["predictors"]][[i]]$parameters)){
        data[[paste0("prior_",paste0("beta_", names(priors$predictors)[i]),"_",par)]] <-  priors[["predictors"]][[i]]$parameters[[par]]
      }
    }
  }

  return(data)
}
.fit_inits             <- function(priors, distribution, chains, seed){

  if(is.null(seed)){
    seed <- sample(666666, 1)
  }
  inits <- vector(mode = "list", chains)
  set.seed(seed)

  for(i in 1:chains){

    temp_init <- list()

    temp_init <- c(temp_init, .fit_inits_mu_tau(priors[["intercept"]], "intercept"))
    if(.has_aux(distribution)){
      temp_init <- c(temp_init, .fit_inits_mu_tau(priors[["aux"]], "aux"))
    }
    for(j in seq_along(priors[["predictors"]])){
      temp_init <- c(temp_init, .fit_inits_mu_tau(
        priors[["predictors"]][[j]],
        paste0("beta_", names(priors[["predictors"]])[j])))
    }

    temp_init[[".RNG.seed"]] <- seed + i
    temp_init[[".RNG.name"]] <- "base::Super-Duper"

    inits[[i]] <- temp_init
  }

  return(inits)
}
.to_monitor            <- function(distribution, predictors){
  c("intercept", if(.has_aux(distribution)) "aux", if(!is.null(predictors)) paste0("beta_", predictors))
}
.fit_model_RoBSA_wrap  <- function(model_syntax, fit_data, fit_inits, monitor_variables, control){
  if(control$silent){
    fit <- callr::r(
      .fit_model_RoBSA,
      args = list(
        model_syntax      = model_syntax,
        fit_data          = fit_data,
        fit_inits         = fit_inits,
        monitor_variables = monitor_variables,
        control           = control
      )
    )
  }else{
    fit <- .fit_model_RoBSA(
      model_syntax      = model_syntax,
      fit_data          = fit_data,
      fit_inits         = fit_inits,
      monitor_variables = monitor_variables,
      control           = control
    )
  }
  return(fit)
}
.fit_model_RoBSA       <- function(model_syntax, fit_data, fit_inits, monitor_variables, control){

  model_call <- list(
    model           = model_syntax,
    data            = fit_data,
    inits           = fit_inits,
    monitor         = monitor_variables,
    n.chains        = control$chains,
    startburnin     = control$burnin,
    startsample     = control$iter,
    adapt           = control$adapt,
    thin            = control$thin,
    raftery.options = if(control$autofit)  list(r = control$max_error) else FALSE,
    psrf.target     = if(control$autofit)  control$max_rhat else Inf,
    max.time        = if(control$autofit)  control$max_time else Inf,
    method          = if(control$parallel) "rjparallel"     else "rjags",
    summarise       = FALSE
  )

  if(control$parallel){
    # the cluster needs to be created manually, because windows don't share the RoBSA JAGS module with the cluster by default
    cl <- parallel::makePSOCKcluster(if(control$chains > control$cores) control$cores else control$chains)
    parallel::clusterCall(cl, function(x) requireNamespace("RoBSA"))
    model_call$cl <- cl
  }else{
    # requires namespace in case that the fit is estimated in a separate R process (for the silent mode)
    requireNamespace("RoBSA")
  }

  if(!is.null(control$seed))set.seed(control$seed)
  fit <- tryCatch(do.call(runjags::autorun.jags, model_call), error = function(e)e)

  if(control$parallel){
    parallel::stopCluster(cl)
  }

  return(fit)
}

### marglik syntax
.marglik_prepare_data  <- function(fit, priors, distribution, data){

  # posterior samples
  samples <- suppressWarnings(coda::as.mcmc(fit)) # gives a warning about merging chains

  ### get parameteres based on the model type
  pars <- NULL
  lb   <- NULL
  ub   <- NULL

  # intercept in auxiliary parameters
  for(par in c("intercept", if(.has_aux(distribution)) "aux")){

    if(priors[[par]]$distribution != "point"){
      # add parameters and truncation
      if(priors[[par]]$distribution == "invgamma"){
        pars <- c(pars, paste0("inv_", par))
        lb   <- c(lb, priors[[par]]$truncation$upper^-1)
        ub   <- c(ub, priors[[par]]$truncation$lower^-1)
      }else{
        pars <- c(pars, par)
        lb   <- c(lb, priors[[par]]$truncation$lower)
        ub   <- c(ub, priors[[par]]$truncation$upper)
      }
    }

  }

  # predictors
  for(i in seq_along(priors[["predictors"]])){

    if(priors[["predictors"]][[i]]$distribution != "point"){
      # add parameters and truncation
      if(priors[["predictors"]][[i]]$distribution == "invgamma"){
        pars <- c(pars, paste0("inv_beta_", names(priors[["predictors"]])[i]))
        lb   <- c(lb, priors[["predictors"]][[i]]$truncation$upper^-1)
        ub   <- c(ub, priors[["predictors"]][[i]]$truncation$lower^-1)
      }else{
        pars <- c(pars, paste0("beta_", names(priors[["predictors"]])[i]))
        lb   <- c(lb, priors[["predictors"]][[i]]$truncation$lower)
        ub   <- c(ub, priors[["predictors"]][[i]]$truncation$upper)
      }
    }
  }


  # keep only the samples needed for the bridge sampling function
  samples <- samples[,pars,drop = FALSE]

  # name the bounds
  names(ub) <- pars
  names(lb) <- pars

  out <- list(
    samples = samples,
    lb      = lb,
    ub      = ub
  )
  return(out)
}
.marglik_function      <- function(samples.row, data, priors, distribution){

   # intercept
  if(priors[["intercept"]]$distribution != "point"){
    if(priors[["intercept"]]$distribution == "invgamma"){
      inv_intercept <- samples.row[[ "inv_intercept" ]]
      intercept     <- 1/inv_intercept
    }else{
      intercept     <- samples.row[[ "intercept" ]]
    }
  }else{
    intercept <- priors[["intercept"]]$parameters$location
  }

  # auxiliary
  if(.has_aux(distribution)){
    if(priors[["aux"]]$distribution != "point"){
      if(priors[["aux"]]$distribution == "invgamma"){
        inv_aux <- samples.row[[ "inv_aux" ]]
        aux     <- 1/inv_aux
      }else{
        aux     <- samples.row[[ "aux" ]]
      }
    }else{
      aux <- priors[["aux"]]$parameters$location
    }
  }else{
    aux = NULL
  }

  # predictors
  beta <- list()
  for(i in seq_along(priors[["predictors"]])){
    temp_beta <- NULL

    if(priors[["predictors"]][[i]]$distribution != "point"){
      if(priors[["predictors"]][[i]]$distribution == "invgamma"){
        temp_beta <- samples.row[[ paste0("inv_beta_", names(priors[["predictors"]])[i]) ]]
        temp_beta <- 1/temp_beta
      }else{
        temp_beta <- samples.row[[ paste0( "beta_", names(priors[["predictors"]])[i]) ]]
      }
    }else{
      temp_beta <- priors[["predictors"]][[i]]$parameters$location
    }

    beta[[names(priors[["predictors"]])[i]]] <- temp_beta
  }


  # compute the linear predictor
  mu <- list()

  # TODO: add "lcent", "icent", "delay" - can be replaced by attribute type in data
  for(type in attr(data ,"type")){
    mu[[type]] <- intercept
    for(i in seq_along(priors[["predictors"]])){
      mu[[type]] <- mu[[type]] + beta[[names(priors[["predictors"]])[i]]] * data[[paste0("x_", type, "_", names(priors[["predictors"]])[i])]]
    }
  }


  ### compute the marginal log_likelihood
  log_lik <- 0

  # intercept
  log_lik <- log_lik + .marglik_distribution(intercept, "intercept", priors[["intercept"]]$distribution, data, priors[["intercept"]]$truncation)

  # auxiliary parameter
  if(.has_aux(distribution)){
    log_lik <- log_lik + .marglik_distribution(aux, "aux", priors[["aux"]]$distribution, data, priors[["aux"]]$truncation)
  }

  # predictors
  for(i in seq_along(priors[["predictors"]])){
    log_lik <- log_lik + .marglik_distribution(
      beta[[names(priors[["predictors"]])[i]]],
      paste0( "beta_", names(priors[["predictors"]])[i]),
      priors[["predictors"]][[i]]$distribution,
      data,
      priors[["predictors"]][[i]]$truncation)
  }

  # data
  for(type in attr(data ,"type")){
    log_lik <- log_lik + sum(.marglik_survival(data[[paste0("t_", type)]], mu[[type]], aux, distribution, type))
  }

  return(log_lik)
}
.marglik_survival      <- function(x, mu, aux, distribution, type){

  # TODO: add other types of censoring
  survival_likelihood <- eval(parse(text = paste0(gsub("-", "_", distribution), "_log_", switch(
    type,
    "event" = "density",
    "rcent" = "survival",
    "lcent" = NULL,
    "icent" = NULL,
    "delay" = NULL
  ))))
  args <- list(
    t   = x,
    eta = mu
  )
  if(.has_aux(distribution)){
    args <- c(args, list(aux))
  }

  log_lik <- do.call(survival_likelihood, args)

  return(log_lik)
}

### default priors settings
get_default_prior_beta_null <- function(){
  prior("spike", list(0))
}
get_default_prior_beta_alt  <- function(){
  prior("normal", list(mean = 0, sd = 1))
}
get_default_prior_intercept <- function(){
  list(
    "exp-aft"     = prior("normal", list(mean = 0, sd = 5)),
    "weibull-aft" = prior("normal", list(mean = 0, sd = 5)),
    "lnorm-aft"   = prior("normal", list(mean = 0, sd = 5)),
    "llogis-aft"  = prior("normal", list(mean = 0, sd = 5)),
    "gamma-aft"   = prior("normal", list(mean = 0, sd = 5))
  )
}
get_default_prior_aux       <- function(){
  list(
    "exp-aft"     = NULL,
    "weibull-aft" = prior("normal", list(mean = 0, sd = 1), list(0, Inf)),
    "lnorm-aft"   = prior("normal", list(mean = 0, sd = 1), list(0, Inf)),
    "llogis-aft"  = prior("normal", list(mean = 0, sd = 1), list(0, Inf)),
    "gamma-aft"   = prior("normal", list(mean = 0, sd = 1), list(0, Inf))
  )
}

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
    "weibull-aft" = exp,
    "lnorm-aft"   = function(x)x,
    "llogis-aft"  = exp,
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

# from RoBSA (with removal of some arguments)
.set_control       <- function(control, chains, iter, burnin, thin, seed, parallel){

  # set the control list
  if(is.null(control)){
    control$autofit         <- FALSE
    control$adapt           <- 1000
    control$bridge_max_iter <- 10000

    control$allow_max_error <- NULL
    control$allow_max_rhat  <- NULL
    control$allow_min_ESS   <- NULL
    control$allow_inc_theta <- FALSE
    control$balance_prob    <- TRUE

    control$silent          <- FALSE

    if(parallel){
      control$cores         <- parallel::detectCores() - 1
    }else{
      control$cores         <- 1
    }

  }else{
    if(is.null(control[["max_error"]])){
      control$max_error       <- .01
    }
    if(is.null(control[["max_rhat"]])){
      control$max_rhat        <- 1.05
    }
    if(is.null(control[["max_time"]])){
      control$max_time        <- Inf
    }
    if(is.null(control[["autofit"]])){
      control$autofit         <- FALSE
    }
    if(is.null(control[["adapt"]])){
      control$adapt           <- 1000
    }
    if(is.null(control[["bridge_max_iter"]])){
      control$bridge_max_iter <- 10000
    }
    if(is.null(control[["allow_max_error"]])){
      control$allow_max_error <- NULL
    }
    if(is.null(control[["allow_max_rhat"]])){
      control$allow_max_rhat  <- NULL
    }
    if(is.null(control[["allow_min_ESS"]])){
      control$allow_min_ESS   <- NULL
    }
    if(is.null(control[["allow_inc_theta"]])){
      control$allow_inc_theta <- FALSE
    }
    if(is.null(control[["balance_prob"]])){
      control$balance_prob    <- TRUE
    }
    if(is.null(control[["silent"]])){
      control$silent          <- FALSE
    }
    if(is.null(control[["cores"]])){
      if(parallel){
        control$cores         <- parallel::detectCores() - 1
      }else{
        control$cores         <- 1
      }
    }

  }

  if(control[["cores"]] > 1){
    parallel <- TRUE
  }

  # add the main MCMC settings
  control$chains    <- chains
  control$iter      <- iter
  control$burnin    <- burnin
  control$thin      <- thin
  control$seed      <- seed
  control$parallel  <- parallel

  .check_control(control)
  return(control)
}
.update_control    <- function(control, control_new, chains, iter, burnin, thin, seed, parallel){

  if(!is.null(control_new)){
    for(n in names(control_new)){
      control[[n]] <- control_new[[n]]
    }
  }

  if(!is.null(chains))   control[["chains"]]   <- chains
  if(!is.null(iter))     control[["iter"]]     <- iter
  if(!is.null(burnin))   control[["burnin"]]   <- burnin
  if(!is.null(thin))     control[["thin"]]     <- thin
  if(!is.null(parallel)) control[["parallel"]] <- parallel
  if(!is.null(seed))     control[["seed"]]     <- seed

  # stop if there is not enough samples planned for autojags package
  .check_control(control)

  return(control)
}
.check_control     <- function(control){
  # check whether only known controls were supplied
  known_controls <- c("chains", "iter", "burnin" , "adapt", "thin" ,"autofit", "max_error", "max_rhat", "max_time", "bridge_max_iter", "allow_max_error", "allow_max_rhat", "allow_min_ESS", "allow_inc_theta", "balance_prob", "silent", "progress_start", "progress_tick", "cores", "seed", "parallel")
  if(any(!names(control) %in% known_controls))stop(paste0("The following control settings were not recognize: ", paste(names(control[!names(control) %in% known_controls]), collapse = ", ")))

  # check whether essential controls were supplied
  if(is.null(control[["chains"]])) stop("Number of chains must be defined.")
  if(is.null(control[["iter"]]))   stop("Number of iterations must be set.")
  if(is.null(control[["burnin"]])) stop("Number of burnin samples must be set.")
  if(is.null(control[["adapt"]]))  stop("Number of adaptation samples must be set.")
  if(is.null(control[["thin"]]))   stop("Thinning of the posterior samples must be set.")

  if(!is.numeric(control[["chains"]]) | !control[["chains"]] >= 1)  stop("At least one chains must be set.")
  if(!is.numeric(control[["iter"]])   | !control[["iter"]] >= 1)    stop("Number of iterations must be a positive number.")
  if(!is.numeric(control[["burnin"]]) | !control[["burnin"]] >= 1)  stop("Number of burnin samples must be a positive number.")
  if(!is.numeric(control[["adapt"]])  | !control[["adapt"]] >= 1)   stop("Number of adaptation samples must be a positive number.")
  if(!is.numeric(control[["thin"]])   | !control[["thin"]] >= 1)    stop("Thinning of the posterior samples must be a positive number.")
  if(!is.logical(control[["parallel"]]))                            stop("The usage of parallel evaluation must be a logical argument.")
  if(!is.numeric(control[["cores"]])  | !control[["cores"]] >= 1)   stop("Number of cores must be a positive number.")
  if(!is.numeric(control[["seed"]])   & !is.null(control[["seed"]]))stop("Seed must be a numeric argument.")

  # stop if there is not enough samples planned for autojags package
  if(control[["iter"]]/control[["thin"]] < 4000)stop("At least 4000 iterations after thinning is required to compute the Raftery and Lewis's diagnostic.")
  if(!control[["chains"]] >= 2)stop("The number of chains must be at least 2 so that convergence can be assessed.")

  # check convergence criteria
  if(control$autofit)if(control[["max_error"]] >= 1 | control[["max_error"]] <= 0)stop("The target maximum MCMC error must be within 0 and 1.")
  if(control$autofit)if(control[["max_rhat"]] <= 1)stop("The target maximum R-hat must be higher than 1.")
  if(!is.null(control[["allow_max_error"]])) if(control[["allow_max_error"]] >= 1 | control[["allow_max_error"]] <= 0)stop("The maximum allowed MCMC error must be within 0 and 1.")
  if(!is.null(control[["allow_max_rhat"]]))  if(control[["allow_max_rhat"]] <= 1) stop("The maximum allowed R-hat must be higher than 1.")
  if(!is.null(control[["allow_min_ESS"]]))   if(control[["allow_min_ESS"]] <= 0)  stop("The minimum allowed ESS must be higher than 0.")

  if(control[["parallel"]]){
    if(!try(requireNamespace("parallel")))stop("parallel package needs to be installed for parallel processing. Run 'install.packages('parallel')'")
  }
  # now taken care of by the evaluation outside of R
  # runjags::runjags.options(silent.jags = control$silent, silent.runjags = control$silent)
}
.runjags.summary   <- function(fit, distribution){

  # returns quantile intervals instead of HPD and produce output on all scales
  invisible(utils::capture.output(summary_fit <- summary(fit, silent.jags = T)))
  model_samples <- suppressWarnings(coda::as.mcmc(fit))

  for(par in rownames(summary_fit)){
    summary_fit[par, "Lower95"] <- stats::quantile(model_samples[,par], .025)
    summary_fit[par, "Upper95"] <- stats::quantile(model_samples[,par], .975)
  }

  return(summary_fit)
}
.fitting_priority  <- function(models){

  predictors    <- sapply(models, function(m)length(m$predictors))
  distributions <- sapply(models, function(m)switch(
    as.character(m$distribution),
    "exp-aft"     = 0,
    "weibull-aft" = 1.5,
    "lnorm-aft"   = 1,
    "llogis-aft"  = 1,
    "gamma-aft"   = 1.5
  ))

  fitting_difficulty <- predictors + distributions

  return(order(fitting_difficulty, decreasing = TRUE))
}
.get_converged_models       <- function(object){

  converged <- NULL

  # basic convergence checks
  for(i in 1:length(object$models)){
    if(any(class(object$models[[i]]$fit) %in% c("simpleError", "error")) | is.infinite(object$models[[i]]$marg_lik$logml) | is.na(object$models[[i]]$marg_lik$logml)){
      converged <- c(converged, FALSE)
    }else{
      converged <- c(converged, TRUE)
    }
  }

  object$models <- object$models[converged]

  # remove models with unsatisfactory performance
  if(!is.null(object$control$allow_max_error) |!is.null(object$control$allow_max_rhat) | !is.null(object$control$allow_min_ESS)){
    diagnostics_summary <- summary.RoBSA(object, type = "models", diagnostics = TRUE, include_theta = object$control$allow_inc_theta)$diagnostics

    # deal with NAs for null models
    diagnostics_summary$"max(MCMC error)"[is.na(diagnostics_summary$"max(MCMC error)")] <- 0
    diagnostics_summary$"max(Rhat)"[is.na(diagnostics_summary$"max(Rhat)")]             <- 0
    diagnostics_summary$"min(ESS)"[is.na(diagnostics_summary$"min(ESS)")]               <- Inf


    if(!is.null(object$control$allow_max_error)){
      converged <- converged & (diagnostics_summary$"max(MCMC error)" < object$control$allow_max_error)
    }
    if(!is.null(object$control$allow_max_Rhat)){
      converged <- converged & diagnostics_summary$"max(Rhat)" < object$control$allow_max_rhat
    }
    if(!is.null(object$control$allow_min_ESS)){
      converged <- converged & diagnostics_summary$"min(ESS)"  > object$control$allow_min_ESS
    }
  }

  return(converged)
}
.model_refit_warnings       <- function(metadata){

  new_warn <- NULL

  # extract meta-data with fit-refit information
  refit_info <- t(sapply(metadata, function(x){
    if(is.null(x$refit_info)){
      return(c(x$i, NA))
    }else{
      return(c(x$i, x$refit_info))
    }
  }))

  marglik_info <- t(sapply(metadata, function(x){
    if(is.null(x$marg_lik)){
      return(c(x$i, NA))
    }else{
      return(c(x$i, x$marg_lik))
    }
  }))

  if(is.null(dim(refit_info)))  refit_info   <- matrix(refit_info,   ncol = 2)
  if(is.null(dim(marglik_info)))marglik_info <- matrix(marglik_info, ncol = 2)


  if(length(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Initial fit of %1$s %2$s failed due to incompatible starting values (most likely due to an outlier in the data and limited precision of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data.",
      ifelse(length(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1]) == 1, "model", "models"),
      paste(refit_info[refit_info[, 2] == "empirical init" & !is.na(refit_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Initial fit of %1$s %2$s failed due to incompatible starting values (most likely due to an outlier in the data and limited precision of t-distribution). Starting values for the mean parameter were therefore set to the mean of supplied data and the model was refitted using boost likelihood function.",
      ifelse(length(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1]) == 1, "model", "models"),
      paste(refit_info[refit_info[, 2] == "refit with boost" & !is.na(refit_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(refit_info[!refit_info[, 2] %in% c("empirical init", "refit with boost") & !is.na(refit_info[,2]), 1]) > 0){
    refit_info_messages_i <- refit_info[refit_info[, 2] != "not enough iterations" & !is.na(refit_info[,2]), 1]
    refit_info_messages   <- refit_info[refit_info[, 2] != "not enough iterations" & !is.na(refit_info[,2]), 2]

    for(i in 1:length(refit_info_messages_i)){
      new_warn <- c(new_warn, paste0("Model ", refit_info_messages_i[i]," failed with the following error: ", refit_info_messages[i]))
    }
  }


  if(length(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1]) > 0){
    new_warn <- c(new_warn, sprintf(
      "Marginal likelihood computation of %1$s %2$s couldn't be completed within the specified number of iterations.",
      ifelse(length(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1]) == 1, "model", "models"),
      paste(marglik_info[marglik_info[, 2] == "not enough iterations" & !is.na(marglik_info[,2]), 1], collapse = ", ")
    ))
  }

  if(length(marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 1]) > 0){
    marglik_info_messages_i <- marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 1]
    marglik_info_messages   <- marglik_info[marglik_info[, 2] != "not enough iterations" & !is.na(marglik_info[,2]), 2]

    for(i in 1:length(marglik_info_messages_i)){
      new_warn <- c(new_warn, paste0("Marginal likelihood computation of model ", marglik_info_messages_i[i]," failed with the following error: ", marglik_info_messages[i]))
    }
  }

  return(new_warn)
}
.model_convergence_warnings <- function(object){

  new_warn <- NULL

  # used set values if specified by the user
  threshold_error <- ifelse(is.null(object$control$allow_max_error), Inf, object$control$allow_max_error)
  threshold_rhat  <- ifelse(is.null(object$control$allow_max_rhat), 1.05, object$control$allow_max_rhat)
  threshold_ESS   <- ifelse(is.null(object$control$allow_max_error), 100, object$control$allow_min_ESS)

  # get the diagnostics summary
  diagnostics_summary <- summary.RoBSA(object, type = "models", diagnostics = TRUE, include_theta = object$control$allow_inc_theta)$diagnostics

  # deal with NAs for null models
  diagnostics_summary$"max(MCMC error)"[is.na(diagnostics_summary$"max(MCMC error)")] <- 0
  diagnostics_summary$"max(Rhat)"[is.na(diagnostics_summary$"max(Rhat)")]             <- 0
  diagnostics_summary$"min(ESS)"[is.na(diagnostics_summary$"min(ESS)")]               <- Inf

  # find the problematic models
  warning_error <- rownames(diagnostics_summary)[diagnostics_summary$"max(MCMC error)" > threshold_error]
  warning_rhat  <- rownames(diagnostics_summary)[diagnostics_summary$"max(Rhat)"       > threshold_rhat]
  warning_ESS   <- rownames(diagnostics_summary)[diagnostics_summary$"min(ESS)"        < threshold_ESS]

  # add warnings messages
  if(length(warning_error) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with MCMC error larger than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_error) == 1, "Model", "Models"),
      paste(warning_error, collapse = ", "),
      threshold_error
    ))
  }

  if(length(warning_rhat) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with R-hat larger than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_rhat) == 1, "Model", "Models"),
      paste(warning_rhat, collapse = ", "),
      threshold_rhat
    ))
  }

  if(length(warning_ESS) > 0){
    new_warn <- c(new_warn, sprintf(
      "%1$s %2$s had at least one parameter with ESS lower than %3$s. We advice checking the MCMC diagnostics before drawing inference from the models or ensemble.",
      ifelse(length(warning_ESS) == 1, "Model", "Models"),
      paste(warning_ESS, collapse = ", "),
      threshold_ESS
    ))
  }

  return(new_warn)
}
.inclusion_BF               <- function(prior_weights, posterior_weights, conditional_models){
  (sum(posterior_weights[conditional_models])/sum(posterior_weights[!conditional_models]))  /
    (sum(prior_weights[conditional_models])/sum(prior_weights[!conditional_models]))
}
.model_inference            <- function(object, n_samples = 10000){

  models    <- object$models
  data      <- object$data
  converged <- object$add_info$converged
  seed      <- object$control$seed

  predictors    <- attr(object$data, "predictors")
  distributions <- as.character(sapply(1:length(object$models), function(i)object$models[[i]]$distribution))

  # extract marginal likelihoods
  marg_liks <- sapply(models, function(x)x$marg_lik$logml)

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
  weights_all        <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_all)
  weights_predictors <- list()
  for(i in seq_along(predictors)){
    if(any(mm_predictors[[predictors[i]]]) & all(!is.nan(prior_weights_predictors[[predictors[i]]]))){
      weights_predictors[[predictors[i]]] <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_predictors[[predictors[i]]])
    }
  }
  weights_distributions <- list()
  for(i in seq_along(unique(distributions))){
    if(any(mm_distributions[[unique(distributions)[i]]]) & all(!is.nan(prior_weights_distributions[[unique(distributions)[i]]]))){
      weights_distributions[[unique(distributions)[i]]] <- bridgesampling::post_prob(marg_liks, prior_prob = prior_weights_distributions[[unique(distributions)[i]]])
    }
  }

  ### update the distributions and predictors trackers only to those present
  predictors    <- names(weights_predictors)
  distributions <- names(weights_distributions)

  ### compute inclusion BFs
  BF_predictors    <- list()
  for(i in seq_along(predictors)){
    BF_predictors[[predictors[i]]] <- .inclusion_BF(prior_weights_all, weights_all, mm_predictors[[predictors[i]]])
  }
  BF_distributions <- list()
  for(i in seq_along(unique(distributions))){
    BF_distributions[[unique(distributions)[i]]] <- .inclusion_BF(prior_weights_all, weights_all, mm_distributions[[unique(distributions)[i]]])
  }


  ### sample and mix the individual posteriors
  if(!is.null(seed))set.seed(seed)
  samples <- list()
  samples$predictors <- list()
  for(i in seq_along(predictors)){
    samples$predictors$averaged[[predictors[i]]]    <- .mix_samples(models, weights_all, converged, predictors[i], n_samples, seed)
  }
  for(i in seq_along(predictors)){
    samples$predictors$conditional[[predictors[i]]] <- .mix_samples(models, weights_predictors[[predictors[i]]], converged, predictors[i], n_samples, seed)
  }
  for(i in seq_along(unique(distributions))){
    samples$intercept[[unique(distributions)[i]]] <- .mix_samples(models, weights_distributions[[unique(distributions)[i]]], converged, "intercept", n_samples, seed)
  }
  for(i in seq_along(unique(distributions))){
    if(.has_aux(unique(distributions)[i])){
      samples$aux[[unique(distributions)[i]]] <- .mix_samples(models, weights_distributions[[unique(distributions)[i]]], converged, "aux", n_samples, seed)
    }
  }


  ### edit names
  names(weights_all)   <- names(models)


  ### compute prior and posterior probabilities for model types
  prior_prob_all            <- prior_weights_all
  prior_prob_distributions  <- sapply(unique(distributions), function(distribution)sum(prior_weights_all[mm_distributions[[distribution]]]))
  prior_prob_predictors     <- sapply(predictors,            function(predictor)   sum(prior_weights_all[mm_predictors[[predictor]]]))

  posterior_prob_all            <- weights_all
  posterior_prob_distributions  <- sapply(unique(distributions), function(distribution)sum(weights_all[mm_distributions[[distribution]]]))
  posterior_prob_predictors     <- sapply(predictors,            function(predictor)   sum(weights_all[mm_predictors[[predictor]]]))


  # return the results
  output <- list(
    samples        = samples,
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
.mix_samples                <- function(models, weights, converged, parameter, n_samples, seed){

  is.predictor <- !parameter %in% c("intercept", "aux")

  if(!is.null(seed)) set.seed(seed) else set.seed(1)
  samples <- NULL

  for(i in c(1:length(models))[converged]){

    model_samples <- suppressWarnings(coda::as.mcmc(models[[i]]$fit))
    # deal with the possibility of intercept only exponential model
    if(!is.matrix(model_samples)){
      model_samples <- matrix(model_samples, ncol = 1)
      colnames(model_samples) <- "intercept"
    }


    # creates indexing - the set seet at the beggining makes sure that the samples for different predictors are correlated
    ind <- sample(nrow(model_samples), round(n_samples * weights[i]), replace = TRUE)
    if(length(ind) == 0)next

    if(is.predictor){
      samples <- c(samples,
                   if(models[[i]]$priors[["predictors"]][[parameter]]$distribution == "point"){
                     rep(models[[i]]$priors[["predictors"]][[parameter]]$parameters$location, round(n_samples * weights[i]))
                   }else{
                     model_samples[ind, paste0("beta_", parameter)]
                   })
    }else{
      samples <- c(samples,
                   if(models[[i]]$priors[[parameter]]$distribution == "point"){
                     rep(models[[i]]$priors[[parameter]]$parameters$location, round(n_samples * weights[i]))
                   }else{
                     model_samples[ind, parameter]
                   })
    }

  }


  return(samples)
}
.compute_coeficients        <- function(RoBSA){
  sapply(RoBSA$samples$predictors$averaged, mean)
}
.BF_format                  <- function(BF, BF01 = FALSE, logBF = FALSE){
  BF[is.nan(BF)] <- NA
  if(BF01){
    BF <- 1/BF
  }else{
    BF <- BF
  }
  if(logBF){
    BF <- log(BF)
  }
  return(BF)
}
