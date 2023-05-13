predict3.bart <- function(object, x.layers, quantiles = NULL) {
  var_names <- attr(object$fit$data@x, "term.labels")
  complete_inds <- which(complete.cases(x.layers[ , var_names]))
  pred_all <- stats::predict(object, x.layers, na.action = na.omit())  # smaller nrow than x.layers if there are NAs
  pred_mean <- colMeans(pred_all)
  if (is.null(quantiles)) {
    out <- NA_real_
    out[complete_inds] <- pred_mean
    return(out)
  } else {
    pred_quants <- apply(pred_all, 2, quantile, probs = quantiles)
    pred_quants <- as.data.frame(t(pred_quants))
    # out <- data.frame(pred = pred_mean, lower = pred_quants[ , 1], upper = pred_quants[ , 2])
    out <- data.frame(matrix(data = NA_real_, nrow = nrow(x.layers), ncol = 3))
    colnames(out) <- c("pred", "lower", "upper")
    out$pred[complete_inds] <- pred_mean
    out$lower[complete_inds] <- pred_quants[ , 1]
    out$upper[complete_inds] <- pred_quants[ , 2]
    return(out)
  }
}
