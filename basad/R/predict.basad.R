### predict.basad.R
### The default predict method for the basad
###
###
### Author: Qingyan Xiang

predict.basad <- function(object, newdata = NULL, ...) {
  #### Check the testData
  if (missing(newdata)) {
    newdata <- object$x
  }

  if (any(is.na(newdata))) {
    stop("NA not permitted in x")
  }

  if (is.null(colnames(newdata))) {
    if ((ncol(newdata) + 1) != ncol(object$x)) stop("test data dimension does not match training data, variable names are not supplied...")
  }

  newdata <- cbind(rep(1, nrow(newdata)), newdata)


  newB <- numeric(ncol(object$x))
  newB[object$model.index] <- object$est.B[object$model.index]
  y <- newdata %*% newB



  return(y)
}
