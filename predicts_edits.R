# functions from the 'predicts' R package by R. Hijmans (https://github.com/rspatial/predicts)
# edited by A.M. Barbosa (https://github.com/AMBarbosa) to accommodate gbm3 and randomForest models
# PR submitted to 'predicts' repo (https://github.com/rspatial/predicts/pull/21)


RMSE <- function(obs, prd, na.rm=FALSE) {
  obs <- as.numeric(obs)
  prd <- as.numeric(prd)
  sqrt(mean((obs - prd)^2, na.rm=na.rm))
}


varImportance <- function(model, data, vars=colnames(data), n=10, ...) {
  RMSE <- matrix(nrow=n, ncol=length(vars))
  colnames(RMSE) <- vars

  if (missing(data)) {
    data <- .get_model_data(model)
    if (is.null(data)) {
      stop("data argument cannot be missing when using this model type")
    }
  }

  P <- predict(model, data, ...)
  for (i in 1:length(vars)) {
    rd <- data
    v <- vars[i]
    for (j in 1:n) {
      rd[[v]] <- sample(rd[[v]])
      p <- predict(model, rd, ...)
      RMSE[j,i] <- RMSE(P, p)
    }
  }
  colMeans(RMSE)
}
