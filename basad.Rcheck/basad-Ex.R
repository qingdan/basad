pkgname <- "basad"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('basad')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("basad")
### * basad

flush(stderr()); flush(stdout())

### Name: basad
### Title: bayesian variable selection with shrinking and diffusing priors
### Aliases: basad
### Keywords: regression

### ** Examples


## Not run: 
##D #------------------------------------------------------
##D Example 1:
##D #------------------------------------------------------
##D 
##D obj <- basad( x = X, y = Y)
##D obj
## End(Not run)

## Not run: 
##D #------------------------------------------------------
##D Example 2: using different priors and slection criteria
##D #------------------------------------------------------
##D 
##D obj <- basad( x = X, y = Y, prior.dist = "t", select.cri = "BIC")
##D obj
##D 
##D 
## End(Not run)






cleanEx()
nameEx("print.basad")
### * print.basad

flush(stderr()); flush(stdout())

### Name: print.basad
### Title: Print Summary Output of Analysis
### Aliases: print.basad
### Keywords: print

### ** Examples

## Not run: 
##D 
##D res <- spikeslab(x = X, y = Y)
##D print(res)
##D 
## End(Not run)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
