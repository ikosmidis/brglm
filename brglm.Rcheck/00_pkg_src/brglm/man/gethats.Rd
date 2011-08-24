\name{gethats}
\alias{gethats}
\title{Calculates the Leverages for a GLM through a C Routine}
\description{
  Calculates the leverages of a GLM through a C routine.  It is intended to
  be used only within \code{\link{brglm.fit}}.
}
\usage{
gethats(nobs, nvars, x.t, XWXinv, ww)
}
\arguments{
  \item{nobs}{The number of observations, i.e. \code{dim(X)[1]}.}
  \item{nvars}{The number of parameters, i.e. \code{dim(X)[1]}, where
  \code{X} is the model matrix, excluding the columns that correspond to
  aliased parameters.}
  \item{x.t}{\code{t(X)}.}
  \item{XWXinv}{The inverse of the Fisher information.}
  \item{ww}{The \sQuote{working} weights.}
}
\value{
  A vector containing the diagonal elements of the hat matrix.
}

\author{Ioannis Kosmidis, \email{ioannis@stats.ucl.ac.uk}}
\seealso{\code{\link{hatvalues}}, \code{\link{brglm.fit}}}
\keyword{regression}
