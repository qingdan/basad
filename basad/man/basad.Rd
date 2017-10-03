\name{basad}
\alias{basad}
\title{Bayesian variable selection with shrinking and diffusing priors}
\description{
    This function performs the Bayesian variable selection procedure with shrinking and diffusing priors via Gibbs sampling. Three different prior options placed on the coefficients are provided: Gaussian, Student's t, Laplace. The posterior estimates of coefficients are returned and the resulted model are selected either by "BIC" criteria or median probability model.
}
\usage{
basad( x = NULL, y = NULL, K = -1, df = 5, nburn = 800, niter = 1500,alternative = FALSE,
    verbose = FALSE, nsplit = 10, prior.dist = "Gauss", select.cri = "median" )
}
\arguments{
\item{x}{The matrix or dataframe of covariates. }
\item{y}{The response variables. }
\item{K}{The initial guess of the numbers of active covariates in the model. This value is related to the prior probability that a covariate is nonzero. If \code{K} is not specified greater than 3, this prior probability will be estimated by a beta prior using Gibbs sampling (see details below).  }
\item{df}{The degree of freedom of t prior when \code{prior.dist == "t"}. }
\item{nburn}{The iteration times of burn-in length (i.e., discarded values). }
\item{niter}{The iteration times of sample length. }
\item{alternative}{Option of whether a fast sampling scheme from Bhattacharya will be used to accelerate the speed of algorithm. However when using block updating( set \code{nsplit} greater than 1) this alternative sampling will not be invoked.  }
\item{verbose}{If TRUE, verbose output is sent to the terminal. }
\item{nsplit}{Numbers of splitting of the block updating. }
\item{prior.dist}{Different prior choices. If \code{prior.dist == "t"}, the algorithm will place t prior for coefficients. If \code{prior.dist == "Laplace"}, the algorithm will place Laplace prior for coefficients. Otherwise it will place the default Gaussian priors. }
\item{select.cri}{Model selection criteria. If \code{select.cri == "median"}, the algorithm will use the median probability model to select the coefficients that are nonzero. If \code{select.cri == "BIC"}, the algorithm will use the BIC criteria to select the coefficients that are nonzero. }
}


\value{

An object of class \code{basad} with the following components:
\item{all.var}{ Summary object for the all variables. }
\item{select.var}{ Summary object for the choosed variables. }
\item{beta.names}{ Variable names for the coefficients. }
\item{verbose}{ Verbose details (used for printing). }
\item{posteriorZ}{A vector of posterior probability of latent variable Z. }
\item{modelIdx}{A vector of index of selected coefficients indicating which coefficients are not zero. }
\item{modelZ}{A binary vector Z indicates whether the coefficient is true in the selected model. }
\item{B}{ Estimated coefficient values from posterior distribution through Gibbs sampling. }
\item{x}{ Standardized x-matrix. }
\item{y}{ Standardized y vector. }
}
\details{
The prior probability \eqn{q_n = P(Z_i = 1)} that a covariate is nonzero can be specified by value K. The K represents a prior belief of the upper bound of the true covariates in the model. When user specify a value of K greater than 3, set \eqn{q_n = c/p_n}, through the calculation( see details in Naveen (2014) ):
\deqn{ \Phi( (K - c)/ \sqrt(c) ) = 1 - \alpha }
we get the prior probability on the models with sizes greater than K will be \eqn{\alpha}, and this \eqn{\alpha} is set to 0.1 in the package.
}
\references{
Narisetty, N. N., & He, X. (2014). Bayesian variable selection with shrinking and diffusing priors. \emph{The Annals of Statistics}, 42(2), 789-817.

Bhattacharya, A., Chakraborty, A., & Mallick, B. K. (2016). Fast sampling with Gaussian scale mixture priors in high-dimensional regression. \emph{Biometrika}, asw042.

Barbieri, M. M., & Berger, J. O. (2004). Optimal predictive model selection. \emph{The annals of statistics}, 32(3), 870-897.
}
\examples{

\dontrun{

#-----------------------------------------------------------
Generate Data: The simulated high dimensional data
#-----------------------------------------------------------

n =  100; p = 499; nz = 5

rho1=0.25;rho2=0.25;rho3=0.25  ### correlations
Bc = c( 0,seq(0.6,3,length.out=nz), array(0, p-nz))

covr1=(1- rho1)*diag(nz) +  array(rho1,c(nz,nz))
covr3=(1- rho3)*diag(p-nz) +  array(rho3,c(p-nz,p-nz))
covr2=array(rho2,c(nz,p-nz))
covr=rbind( cbind(covr1,covr2), cbind(t(covr2),covr3) )

covE = eigen(covr)
covsq = covE$vectors \%*\% diag(sqrt(covE$values)) \%*\% t(covE$vectors)

Xs = matrix( rnorm(n*p), nrow = n); Xn = covsq%*% t(Xs)
X = cbind(array(1, n), t(Xn))
Y = X \%*\% Bc + rnorm(n); X <- X[,2:ncol(X)]


#-----------------------------------------------------------
Example 1: Run the default setting the Guassian prior
#-----------------------------------------------------------

obj <- basad( x = X, y = Y)
print( obj )

#-----------------------------------------------------------
Example 2: using different priors and slection criteria
#-----------------------------------------------------------

obj <- basad( x = X, y = Y, prior.dist = "t", select.cri = "BIC")
print( obj )
}



}
\author{
Qinyan Xiang (\email{qyxiang@illinois.edu})

Naveen Narisetty ( \email{naveen@illinois.edu} )
}
\keyword{regression}

