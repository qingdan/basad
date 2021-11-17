\name{summary.basad}
\alias{summary.basad}
\title{Print summary output of basad}
\description{
    Generate summaries from \code{basad} function.
    This function allows for the choice of selection criterion (median probability model, BIC) to perform a variable selection.}
\usage{
\method{summary}{basad}(object, select.cri = "median", BIC.maxsize = 20, ...)
}
\arguments{
    \item{object}{An object of class \code{basad}.}
    \item{select.cri}{Model selection criterion. If \code{select.cri="median"}, the algorithm will use the median probability model to select the active variables. If \code{ select.cri="BIC"}, the algorithm will use the BIC criterion to select the active variables. }
    \item{BIC.maxsize}{ The maximum number of variables that are chosen to apply BIC criterion based on the ranking of their marginal posterior probabilities. If the dimension is less than the default value of 20, all the variables will be considered to apply BIC. }
    \item{...}{Further arguments passed to or from other methods.}

}
\references{
Narisetty, N. N., & He, X. (2014). Bayesian variable selection with shrinking and diffusing priors. \emph{The Annals of Statistics}, 42(2), 789-817.
}
\author{
    Qingyan Xiang (\email{qyxiang@bu.edu})

    Naveen Narisetty (\email{naveen@illinois.edu})
}
\keyword{summary}
