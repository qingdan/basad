\name{predict.basad}
\alias{predict.basad}
\title{basad prediction}
\description{
    Print summary output from spike and slab analysis.
    Note that this is the default print method for the package.
}
\usage{
\method{predict}{basad}(object, testx, ...)
}
\arguments{
    \item{object}{An object of class \code{basad}.}
    \item{testx}{Data frame or x-matrix containing test data.}
	\item{...}{Further arguments passed to or from other methods.}

}
\details{
 	Computes the predicted value using a test data set.
}
\value{
	A vector of fitted values for basad estimated beta values.
}
\references{
Narisetty, N. N., & He, X. (2014). Bayesian variable selection with shrinking and diffusing priors. \emph{The Annals of Statistics}, 42(2), 789-817.
}
\examples{

}
\author{
Qinyan Xiang (\email{qyxiang16@gmail.com})

Naveen Narisetty ( \email{naveen@illinois.edu} )

}
\keyword{regression}
