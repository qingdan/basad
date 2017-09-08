###
###
###The default print method for the basad
###
###


print.basad <- function(x, ...)
{

#### Check that object is compatible
  if (!inherits(x, "basad"))
     stop("This function only works for objects of class 'basad'")

#### extract summary data
  verboseList = x$verbose
  
#### --------------------------------------------------------------
###	Terminal Output
### --------------------------------------------------------------	


### basad output
    
cat("-----------------------------", "\n")
cat("Sample size                    :", verboseList[[1]], "\n" )
cat("No. predictors                 :", verboseList[[2]], "\n" )
cat("Burn-in periods                :", verboseList[[3]], "\n" )
cat("Sampled periods                :", verboseList[[4]], "\n" )
cat("Alternative Sampling                  :", verboseList[[5]], "\n" )
cat("Model selection criteria       :", verboseList[[6]], "\n")
cat("\n\n")
cat("---> Top variables:\n")
print(round(x$basad.sum, 3))
cat("-----------------------------", "\n")


}
