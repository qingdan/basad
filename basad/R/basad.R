###
###
###Main fucntion of basad
###
###

basad <- function(x = NULL,
                  y = NULL,
                  K = -1,
                  df = 5,
                  nburn = 500,
                  niter = 3000,
                  Fast = TRUE,
                  verbose = FALSE,
                  nsplit = 10,
                  prior.dist = "Gauss",
                  select.cri = "median"){

    #######Check the data
	if (is.null(y) | is.null(x))
		stop("x and or y is missing")
	
	Y <- y
	X <- x
    X <- cbind( rep(1, nrow(X)), X)
    beta.names <- colnames(X)
    if (length(unique(beta.names)) != ncol(X)) {
		colnames(X) <- beta.names <- paste("x.", 1:ncol(X), sep = "")
        }
	colnames(X)[1] <- beta.names[1] <- "intercept"			  
	
    X <- data.matrix(X)
    Y <- data.matrix(Y)
				  
    p = dim(X)[2]-1;  n = dim(X)[1]


    ########calculate the prior probability of q
    
    if( K > 0 ){
        choicep <-function(x){
            return(x - K + qnorm(0.9)*sqrt(x*(1- x/p)))
        }
        cp = uniroot(choicep, c(1,K))$root
    
        pr = cp/p
    }
    else
        pr = -1
    
    pr0 = 0.1
    
    
    ##########calculate the ols estimate of sig
    Bols = sapply(seq(1:(p+1)), function(j) summary(lm(Y~ 0+X[,j]))$coef[,1])
    B0 = Bols
    if(p > n) {m = round(p - (n/2))}
    if(p <=n) m = round(p/2)
    ind = which(abs(B0)> sort(abs(B0))[m+1])
    sighat = sum((summary(lm(Y~0+X[,union(1,ind)]))$residuals)^2)/(n-length(union(1,ind)))

    B0 = rep(0,(p+1))
    Z0 = array(0,(p+1))
    sig = sighat
	
	cat("Algorithms running:",  "\n" )


    #######excute the algorithm
    if( prior.dist == "t"){
        
        s0 = 1/n
        nu = df
        fvalue = dt(sqrt(2.1*log(p+1)), df = nu)
        s1  = max(100*s0, pr0*s0/((1 - pr0)*fvalue));
        
        
		res <- .Call( 'basadFunctionT', X, Y, Z0, B0, sig, pr, n, p, nu, s0, s1, nburn, niter, nsplit, Fast, PACKAGE = 'basad' )
    }
    else if( prior.dist == "Laplace"){
        
        s0 = 1/n
        fvalue = dlaplace(sqrt(2.1*log(p+1)))
        s1  = max(100*s0, pr0*s0/((1 - pr0)*fvalue));
        
	    res <- .Call( 'basadFunctionL', X, Y, Z0, B0, sig, pr, n, p, lambda = 1, s0, s1, nburn, niter, nsplit, Fast, PACKAGE = 'basad' )
    }
    else if( prior.dist == "Gauss"){
    
        s0 = 1/n
        fvalue = dnorm(sqrt(2.1*log(p+1)))
        s1  = max(100*s0, pr0*s0/((1 - pr0)*fvalue));
    
        res <- .Call( 'basadFunctionG', X, Y, Z0, B0, sig, pr, n, p, s0, s1, nburn, niter, nsplit, Fast,  PACKAGE = 'basad')
    }
    else{
       stop("No such prior type, please choose from 'Gauss', 't', 'Laplace'")
    }



	
###---------------------------
###
###SUMMARY
###
###---------------------------
	
	
	



    ZZ <- res$Z[(nburn+1):(nburn+niter), ]
    BB <- res$B[(nburn+1):(nburn+niter), ]
    
    
    #posterior Z and B vector
    Z <- apply(ZZ, 2, mean)
    B <- apply(BB, 2, mean)
	
    modelIdx <- c()
    
    
    ####return the results if selection criteria is BIC
    if( select.cri == "BIC" ){

        Zsort <- sort.int(Z, decreasing = TRUE, index.return = TRUE )
        bic <- numeric(15)
        
        
        for( i in 1:15 ){
            idx <- Zsort$ix[ 1:i ]
            Bbic <- numeric(p + 1)
            Bbic[idx] <- B[idx]
            bic[i] <- n * log( ( t( Y - X %*% Bbic ) %*% ( Y - X %*% Bbic ) ) / n  ) + log(n) * i
        }
        
        minIdx <- which.min(bic)
        modelIdx <-  Zsort$ix[1:minIdx]
        
    }
    else{
        ##median probability model vector
        modelIdx <- which( Z > 0.5 )
        select.cri = "Median Probability Model"
    }
    
    
    returnZ <- numeric(p+1)
    returnZ[modelIdx] <- 1


	
	basad.sum <- B[modelIdx]
	basad.sum <- as.data.frame(  basad.sum )
	rownames( basad.sum ) <- beta.names[modelIdx]
	colnames( basad.sum ) <- "estimated values"
	
    

	
	
	
	
###---------------------------
###
###PRINT
###
###---------------------------


catList <- list(

    c(n),
    c(p+1),
    c(nburn),
    c(niter),
	c(Fast),
    c(select.cri)
    
)

if( verbose ){
cat("-----------------------------", "\n")
cat("Sample size                    :", catList[[1]], "\n" )
cat("No. predictors                 :", catList[[2]], "\n" )
cat("Burn-in periods                :", catList[[3]], "\n" )
cat("Sampled periods                :", catList[[4]], "\n" )
cat("Fast Sampling                  :", catList[[5]], "\n" )
cat("Model selection criteria       :", catList[[6]], "\n" )
cat("\n\n")
cat("---> Top variables:\n")
print(round(basad.sum, 3))
cat("-----------------------------", "\n")
}

###---------------------------
###
###RETURN
###
###---------------------------
    out <- list(
        basad.summary = basad.sum,
		catList = catList,
        n = n,
        p = p+1,
		posteriorZ = Z,
        modelIdx = modelIdx,
        modelZ = returnZ,
		B = B,
		x = X,
		y = Y,
        pr = res$Pr
       )
    
    class(out) <- "basad"
    return( out )
}

