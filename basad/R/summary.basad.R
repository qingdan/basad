### print.basad.R
### The default summary method for the basad
###
###
### Author: Qingyan Xiang


summary.basad <- function(object, select.cri = "median", BIC.maxsize = 20, ...) {
  #### Check that object is compatible
  if (!inherits(object, "basad")) {
    stop("This function only works for objects of class 'basad'")
  }

  X <- object$x
  Y <- object$y
  p <- dim(X)[2] - 1
  n <- dim(X)[1]

  nburn <- object$verbose[[3]]
  niter <- object$verbose[[4]]

  ZZ <- object$allZ[(nburn + 1):(nburn + niter), ]
  BB <- object$allB[(nburn + 1):(nburn + niter), ]


  #### posterior Z and B vector
  Z <- apply(ZZ, 2, mean)
  B <- apply(BB, 2, mean)

  Zsort <- sort.int(Z, decreasing = TRUE, index.return = TRUE)
  modelIdx <- c()

  if (p > BIC.maxsize) {
    BICsize <- BIC.maxsize
  } else {
    BICsize <- p
  }

  if (select.cri == "BIC") {
    print(BICsize)
    bic <- numeric(BICsize)
    for (i in 1:BICsize) {
      idx <- Zsort$ix[1:i]
      Bbic <- numeric(p + 1)
      Bbic[idx] <- B[idx]
      bic[i] <- n * log((t(Y - X %*% Bbic) %*% (Y - X %*% Bbic)) / n) + log(n) * i
    }
    minIdx <- which.min(bic)
    modelIdx <- Zsort$ix[1:minIdx]
  } else {
    #### median probability model vector
    modelIdx <- which(Z > 0.5)
    select.cri <- "Median Probability Model"
  }

  ## rerturn all variables
  basad.all <- round(B[Zsort$ix], 4)
  posterior.all <- round(Z[Zsort$ix], 4)
  basad.all <- data.frame(basad.all, posterior.all)

  rownames(basad.all) <- object$beta.names[Zsort$ix]
  colnames(basad.all) <- c("estimated values", "posterior prob")

  ## return slected variables
  basad.select <- round(B[modelIdx[-1]], 4)
  posterior.select <- round(Z[modelIdx[-1]], 4)
  basad.select <- data.frame(basad.select, posterior.select)
  rownames(basad.select) <- object$beta.names[modelIdx[-1]]
  colnames(basad.select) <- c("estimated values", "posterior prob")


  ## return selected model
  modelZ <- numeric(p + 1)
  modelZ[modelIdx] <- 1
  names(modelZ) <- object$beta.names

  object$select.var <- basad.select
  object$modelSize <- length(modelIdx) - 1
  object$model.index <- modelIdx
  object$modelZ <- modelZ

  #### extract summary data
  verboseList <- object$verbose
  verbostList[[7]] <- select.cri


  #### --------------------------------------------------------------
  ### 	Terminal Output
  ### ---------------------------------------------------------------
  ### basad output
  cat("\nCall: ", deparse(object$call), "\n\n")


  cat("----------------------------------", "\n")
  cat("Sample size                      :", verboseList[[1]], "\n")
  cat("Dimension                        :", verboseList[[2]], "\n")
  cat("Burn-in length                   :", verboseList[[3]], "\n")
  cat("Iteration length                 :", verboseList[[4]], "\n")
  cat("Block updating split sizes       :", verboseList[[5]], "\n")
  cat("Alternative fast sampling        :", verboseList[[6]], "\n")
  cat("Model selection criterion        :", verboseList[[7]], "\n")
  cat("\n\n")
  cat("-----> Selected variables:\n")
  print(round(object$select.var, 4))
  cat("----------------------------------", "\n")
}
