\name{glm.control1}
\alias{glm.control1}
\title{Auxiliary for Controlling BRGLM Fitting}
\description{
  Auxiliary function as user interface for \code{\link{brglm}}
  fitting. Typically only used when calling \code{brglm} or
  \code{brglm.fit}.
}
\usage{
glm.control1(epsilon = 1e-08, maxit = 25, trace = FALSE, ...)
}
\arguments{
  \item{epsilon}{as in \code{\link{glm.control}}.}
  \item{maxit}{as in \code{\link{glm.control}}.}
  \item{trace}{as in \code{\link{glm.control}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The only difference with \code{\link{glm.control}} is that
  \code{glm.control1} supports further arguments to be passed from other
  methods. However, this additional arguments have no effect on the
  resultant list.
}
\author{Ioannis Kosmidis, \email{i.kosmidis@ucl.ac.uk}}
\keyword{iteration}