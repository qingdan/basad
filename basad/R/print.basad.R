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
  catList = x$catList
  
#### --------------------------------------------------------------
###	Terminal Output
### --------------------------------------------------------------	


### basad output
    
cat("-----------------------------", "\n")
cat("Sample size                    :", catList[[1]], "\n" )
cat("No. predictors                 :", catList[[2]], "\n" )
cat("Burn-in periods                :", catList[[3]], "\n" )
cat("Sampled periods                :", catList[[4]], "\n" )
cat("Fast Sampling                  :", catList[[5]], "\n" )
cat("Model selection criteria       :", catList[[6]], "\n")
cat("\n\n")
cat("---> Top variables:\n")
print(round(x$basad.sum, 3))
cat("-----------------------------", "\n")


}
