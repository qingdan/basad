\name{basad}
\alias{basad}
\title{bayesian variable selection with shrinking and diffusing priors}
\description{
    This function performs the bayesian variable selection procedure with shrinking and diffusing priors via Gibbs sampling. Three different prior options placed on the coefficients are provided: guassian, student's t, Laplace. Through "BIC" criteria or median probability model the function returns the variable selection results.
}
\usage{
basad( x = NULL, y = NULL, K = -1, df = 5, nburn = 500, niter = 3000,
    Fast = TRUE, verbose = FALSE, prior.dist = "Gauss", select.cri = "median" )
}
\arguments{
  \item{x}{The matrix or dataframe of covariates.}

\item{y}{The response variables.}

\item{K}{The initial guess of the active covariates that is related to the prior probability of whether a covariate is not zero. If \code{K} is not specified as a positive value, this prior probability will be estimated by a beta prior using Gibbs sampling.}

\item{df}{The degree of freedom of t prior when \code{prior.dist == "t"}. }

\item{nburn}{The iteration times of burning period (i.e., discarded values). }

\item{niter}{The iteration times after burning period.}

\item{Fast}{Option whether a faster sampling scheme from Bhattacharyya will be used to accelerate the algorithm, the default value is \code{TRUE}. }
\item{verbose}{If TRUE, verbose output is sent to the terminal.}
\item{prior.dist}{ Different prior choices, if \code{prior.dist == "t"}, the algorithm will place t prior for coefficients, if \code{prior.dist == "Lap"}, the algorithm will palce Laplace prior for coefficients. Otherwise it will place the default Gaussian priors.  }
\item{select.cri}{ Model selection criteria,if \code{select.cri == "median"}, the algorithm will use the median probability model to select the coefficients that are not zero, if \code{select.cri == "BIC"}, the algorithm will use the BIC criteria to select the coefficients that are not zero.  }

}
\value{

An object of class \code{basad} with the following components:
\item{basad.summary}{Summary object for the choosed variables.}
\item{catList}{Verbose details (used for printing).}
\item{n}{ The number of observations }
\item{p + 1}{ The dimensions of the predictors }
\item{posteriorZ}{ Vector of posterior probability of Z.}
\item{modelIdx}{ A vector of index that which coeffcients are not zero thus selected.}
\item{modelZ}{ Binary vector Z that idicating whether the coffecient is true in the model.}
\item{B}{Coefficients results from Gibbs Sampling.}
\item{x}{Original x-matrix.}
\item{y}{Original y vector.}
}
\details{

}
\references{
Narisetty, N. N., & He, X. (2014). Bayesian variable selection with shrinking and diffusing priors. \emph{The Annals of Statistics}, 42(2), 789-817.

Barbieri, M. M., & Berger, J. O. (2004). Optimal predictive model selection. \emph{The annals of statistics}, 32(3), 870-897.

Bhattacharya, A., Chakraborty, A., & Mallick, B. K. (2016). Fast sampling with Gaussian scale mixture priors in high-dimensional regression. \emph{Biometrika}, asw042.

}
\examples{

\dontrun{
#------------------------------------------------------
Example 1:
#------------------------------------------------------

obj <- basad( x = X, y = Y)
obj
}

\dontrun{
#------------------------------------------------------
Example 2: using different priors and slection criteria
#------------------------------------------------------

obj <- basad( x = X, y = Y, prior.dist = "t", select.cri = "BIC")
obj


}



}
\author{
Qinyan Xiang (\email{qyxiang16@gmail.com})

Naveen Narisetty ( \email{naveen@illinois.edu} )
}
\keyword{regression}

