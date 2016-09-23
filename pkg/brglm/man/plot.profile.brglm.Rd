\name{plot.profile.brglm}
\alias{plot.profile.brglm}
\alias{pairs.profile.brglm}
\title{Plot methods for 'profile.brglm' objects}
\description{
  \code{plot.profile.brglm} plots the objects of class
  \code{"profileModel"} that  are contained in an object of class
  \code{"profile.brglm"}. \code{pairs.profile.brglm} is a diagnostic tool
  that plots pairwise profile traces.
}
\usage{
\method{plot}{profile.brglm}(x, signed = FALSE, interpolate = TRUE,
     n.interpolations = 100, print.grid.points = FALSE, ...)

\method{pairs}{profile.brglm}(x, colours = 2:3, ...)
}
\arguments{
  \item{x}{a \code{"profile.brglm"} object.}
  \item{signed}{as in \code{\link[profileModel]{plot.profileModel}}.}
  \item{interpolate}{as in \code{\link[profileModel]{plot.profileModel}}.}
  \item{n.interpolations}{as in \code{\link[profileModel]{plot.profileModel}}.}
  \item{print.grid.points}{as in \code{\link[profileModel]{plot.profileModel}}.}
  \item{colours}{as in \code{\link[profileModel]{plot.profileModel}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  See Details in \code{\link[profileModel]{plot.profileModel}}.
}
\author{Ioannis Kosmidis, \email{i.kosmidis@ucl.ac.uk}}
\seealso{\code{\link[profileModel]{plot.profileModel}}, \code{\link{profile.brglm}}.}
\examples{
# see example in 'confint.brglm'.
}
\keyword{dplot}
\keyword{hplot}
