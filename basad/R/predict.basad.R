### predict.basad.R
### The default predict method for the basad
###
###
### Author: Qingyan Xiang

predict.basad <- function(object, testx = NULL, ...)
{
    
    
    #### Check the testData
	if (missing(testx))
        stop("testing data missing ...")
    if( any( is.na(testx)) )
        stop("NA not permitted in x")
    
    if (is.null(colnames(testx)))
    {
    if ( (  ncol(testx) + 1 ) != object$p) stop("test data dimension does not match training data, variable names are not supplied...")
    }else if (any(colnames(testx) != colnames( object$x ) ))
    {
        warning("test data variables names does not match training data...")
        varmatch = match( colnames( object$x ), colnames(testx))
        if (any(is.na(varmatch))) stop("test data missing some variables...")
            testx = testx[, varmatch]
    }
    
    
    testx = cbind( rep(1, nrow(testx)), testx )


    newB <- numeric( object$p )
    newB[object$modelIdx] <- object$B[object$modelIdx]
    y = testx %*% newB
    
    return(y)
}
