### print.basad.R
### The default print method for the basad
###
###
### Author: Qingyan Xiang


coef.basad <- function(object, ...) {
  #### extract summary data
  verboseList <- object$verbose

  #### --------------------------------------------------------------
  ### 	Terminal Output
  ### ---------------------------------------------------------------

  print(round(object$all.var, 4))
}
