\name{brglm.control}
\alias{brglm.control}
\title{Auxiliary for Controlling BRGLM Fitting}
\description{
Auxiliary function as user interface for \code{\link{brglm}}
fitting. Typically only used when calling \code{brglm} or \code{brglm.fit}.
}
\usage{
brglm.control(br.epsilon = 1e-08, br.maxit = 100, br.trace=FALSE,
              br.consts = NULL, ...)
}
\arguments{
  \item{br.epsilon}{positive convergence tolerance for the iteration
    described in \code{\link{brglm.fit}}.}
  \item{br.maxit}{integer giving the maximum number of iterations for
  the iteration in \code{\link{brglm.fit}}.}
  \item{br.trace}{logical indicating if output should be prooduced for
    each iteration.}
  \item{br.consts}{a (small) positive constant or a vector of
    such.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  If \code{br.trace=TRUE} then for each iteration the iteration number
  and the current value of the modified scores is
  \code{\link{cat}}'ed. If \code{br.consts} is specified then \code{br.consts}
  is added to the original binomial counts and \code{2*br.consts}. Then
  the model is fitted to the adjusted data to provide starting values
  for the iteration in \code{\link{brglm.fit}}. If \code{br.consts = NULL}
  (default) then \code{\link{brglm.fit}} adjusts the responses and totals by
  "number of parameters"/"number of observations" and twice that, respectively.
}
\value{
  A list with the arguments as components.
}
\references{
  Kosmidis, I. (2007). Bias reduction in exponential family nonlinear
  models. \emph{PhD Thesis}, Department of Statistics, University of
  Warwick.
}
\author{Ioannis Kosmidis, \email{i.kosmidis@ucl.ac.uk}}
\seealso{
  \code{\link{brglm.fit}}, the fitting procedure used by
  \code{\link{brglm}}.
}
\keyword{iteration}
