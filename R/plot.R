#' @title Plots a fitted RoBSA object
#'
#' @param x a fitted RoBSA object.
#' @param ... additional arguments.
#' @export
plot.RoBSA <- function(x, ...){

  # create times for the prediction
  df          <- x$data
  time_range  <- c(0, max(c(df$t_event, df$t_rcent)))
  pred_data   <- data.frame(time = seq(time_range[1], time_range[2], length.out = 100))

  # add mean values for the predictors
  pred_names  <- attr(df, "predictors")
  for(i in seq_along(pred_names)){
    mean_pred   <- mean(unlist(df[grep(pred_names[i], names(df))]))
    pred_data   <- cbind(pred_data, new_pred = rep(mean_pred, nrow(pred_data)))
    colnames(pred_data)[ncol(pred_data)] <- pred_names[i]
  }

  # predict the mean survival
  pred_survival <- predict.RoBSA(x, pred_data, type = "survival")

  # plot the survival
  graphics::plot(x = pred_data$time, y = pred_survival$mean, ylim = c(0,1), "l", xlab = "Time", ylab = "Survival", las = 1, ...)
  graphics::lines(x = pred_data$time, y = pred_survival$uCI, lty = 2, ...)
  graphics::lines(x = pred_data$time, y = pred_survival$lCI, lty = 2, ...)

  return(invisible())
}
