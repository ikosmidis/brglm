\name{profileObjectives-brglm}
\alias{profileObjectives}
\alias{penalizedDeviance}
\alias{modifiedScoreStatistic}
\title{Objectives to be profiled}
\description{
  Objectives that are used in \code{\link{profile.brglm}}
}
\usage{
penalizedDeviance(fm, X, dispersion = 1)

modifiedScoreStatistic(fm, X, dispersion = 1)
}
\arguments{
  \item{fm}{the \bold{restricted} fit.}
  \item{X}{the model matrix of the fit on all parameters.}
  \item{dispersion}{the dispersion parameter.}
}
\details{
  These objectives follow the specifications for objectives in the
  \pkg{profileModel} package and are used from \code{profile.brglm}.

  \code{penalizedDeviance} returns a deviance-like value corresponding
  to a likelihood function penalized by Jeffreys invariant prior. It
  has been used by Heinze & Schemper (2002) and by Bull et. al. (2002)
  for the construction of confidence intervals for the bias-reduced
  estimates in logistic regression. The \code{X} argument is the
  model  matrix of the full (\bold{not} the restricted) fit.
  
  \code{modifiedScoreStatistic} mimics \code{\link[profileModel]{RaoScoreStatistic}}
  in \pkg{profileModel}, but with the ordinary scores replaced with the
  modified scores used for bias reduction. The argument \code{X} has
  the same interpretation as for \code{penalizedDeviance}.
}
\value{
  A scalar.
}
\references{
  Bull, S. B., Lewinger, J. B. and Lee, S. S. F. (2007). Confidence
  intervals for multinomial logistic regression in sparse data.
  \emph{Statistics in Medicine} \bold{26}, 903--918.

  Heinze, G. and Schemper, M. (2002). A solution to the problem of
  separation in logistic regression.
  \emph{Statistics in Medicine} \bold{21}, 2409--2419.
}
\author{Ioannis Kosmidis, \email{I.Kosmidis@warwick.ac.uk}}
\seealso{\code{\link[profileModel]{profileModel}}, \code{\link{profile.brglm}}.}
\keyword{models}