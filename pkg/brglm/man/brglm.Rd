\name{brglm}
\alias{brglm}
\alias{brglm.fit}
\alias{print.brglm}
\alias{summary.brglm}
\alias{print.summary.brglm}
\title{Bias reduction in Binomial-response GLMs}
\description{
Fits binomial-response GLMs using the bias-reduction method developed in
Firth (1993) for the removal of the leading (\eqn{\mathop{\rm
O}(n^{-1})}{O(n^{-1})}) term from the asymptotic expansion of the bias
of the maximum likelihood estimator. Fitting is performed using
pseudo-data representations, as described in Kosmidis (2007, Chapter 5). For
estimation in binomial-response GLMs, the bias-reduction method is an
improvement over traditional maximum likelihood because:
\itemize{
\item the bias-reduced estimator is second-order unbiased and has smaller
      variance than the maximum likelihood estimator and
\item the resultant estimates and their corresponding standard errors
      are \bold{always} finite while the maximum likelihood estimates
      can be infinite (in situations where complete or quasi separation
      occurs).
}
}
\usage{
brglm(formula, family = binomial, data, weights, subset, na.action,
      start = NULL, etastart, mustart, offset,
      control.glm = glm.control1(...), model = TRUE, method = "brglm.fit",
      pl = FALSE, x = FALSE, y = TRUE, contrasts = NULL,
      control.brglm = brglm.control(...), ...)

brglm.fit(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
          mustart = NULL, offset = rep(0, nobs), family = binomial(),
          control = glm.control(), control.brglm = brglm.control(),
          intercept = TRUE, pl = FALSE)
}
\arguments{
  \item{formula}{as in \code{\link{glm}}.}
  \item{family}{as in \code{\link{glm}}. \code{brglm} currently
  supports only the \code{"binomial"} family with links
  \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"cauchit"}.}
  \item{data}{as in \code{\link{glm}}.}
  \item{weights}{as in \code{\link{glm}}.}
  \item{subset}{as in \code{\link{glm}}.}
  \item{na.action}{as in \code{\link{glm}}.}
  \item{start}{as in \code{\link{glm}}.}
  \item{etastart}{as in \code{\link{glm}}.}
  \item{mustart}{as in \code{\link{glm}}.}
  \item{offset}{as in \code{\link{glm}}.}
  \item{control.glm}{\code{control.glm} replaces the \code{control}
  argument in \code{\link{glm}} but essentially does the same job. It
  is a list  of parameters to control \code{\link{glm.fit}}. See the
  documentation of  \code{glm.control1} for details.}
  \item{control}{same as in \code{\link{glm}}. Only available to
    \code{brglm.fit}.}
  \item{intercept}{as in \code{\link{glm}}.}
  \item{model}{as in \code{\link{glm}}.}
  \item{method}{the method to be used for fitting the model. The default
  method is \code{"brglm.fit"}, which uses either the modified-scores
  approach to  estimation or maximum penalized likelihood (see
  the  \code{pl} argument below).  The standard
  \code{\link{glm}} methods \code{"glm.fit"} for maximum likelihood and
  \code{"model.frame"} for returning the model frame without any
  fitting, are also accepted.}
  \item{pl}{a logical value indicating whether the
  model should be fitted using maximum penalized likelihood, where the
  penalization is done using Jeffreys invariant prior, or using the
  bias-reducing modified scores. It is only used when
  \code{method = "brglm.fit"}. The default value is \code{FALSE} (see also the
  Details section).}
  \item{x}{as in \code{\link{glm}}.}
  \item{y}{as in \code{\link{glm}}.}
  \item{contrasts}{as in \code{\link{glm}}.}
  \item{control.brglm}{a list of parameters for controlling the fitting
  process when \code{method = "brglm.fit"}. See documentation of
  \code{\link{brglm.control}} for details.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{brglm.fit} is the workhorse function for fitting the model using
either the bias-reduction method or maximum penalized likelihood. If
\code{method = "glm.fit"}, usual maximum likelihood is used via
\code{\link{glm.fit}}.

The main iteration of \code{brglm.fit} consists of the following
steps:
\enumerate{
\item  Calculate the diagonal components of the hat matrix (see
\code{\link{gethats}} and \code{\link{hatvalues}}).
\item Obtain the pseudo-data representation at the current value of the
parameters (see \code{\link{modifications}} for more information).
\item Fit a local GLM, using \code{\link{glm.fit}} on the pseudo data.
\item Adjust the quadratic weights to agree with the original
binomial totals.
}
Iteration is repeated until either the iteration limit has been reached
or the sum of the absolute values of the modified scores is less than
some specified positive constant (see the \code{br.maxit} and
\code{br.epsilon} arguments in \code{\link{brglm.control}}).

The default value (\code{FALSE}) of \code{pl}, when \code{method = "brglm.fit"},
results in estimates that are free of any \eqn{\mathop{\rm
O}(n^{-1})}{O(n^{-1})} terms in the asymptotic expansion of their bias. When
\code{pl = TRUE} bias-reduction is again achieved but generally not at
such order of magnitude. In the case of logistic regression the value of
\code{pl} is irrelevant since maximum penalized likelihood and the
modified-scores approach coincide for natural exponential families (see
Firth, 1993).

For other language related details see the details section in
\code{\link{glm}}.
}
\value{
  \code{\link{brglm}} returns an object of class \code{"brglm"}. A
  \code{"brglm"} object inherits first from \code{"glm"} and then from
  \code{"lm"} and is a list containing the following components:
  \item{coefficients}{as in \code{\link{glm}}.}
  \item{residuals}{as in \code{\link{glm}}.}
  \item{fitted.values}{as in \code{\link{glm}}.}
  \item{effects}{as in \code{\link{glm}}.}
  \item{R}{as in \code{\link{glm}}.}
  \item{rank}{as in \code{\link{glm}}.}
  \item{qr}{as in \code{\link{glm}}.}
  \item{family}{as in \code{\link{glm}}.}
  \item{linear.predictors}{as in \code{\link{glm}}.}
  \item{deviance}{as in \code{\link{glm}}.}
  \item{aic}{as in \code{\link{glm}} (see Details).}
  \item{null.deviance}{as in \code{\link{glm}}.}
  \item{iter}{as in \code{\link{glm}}.}
  \item{weights}{as in \code{\link{glm}}.}
  \item{prior.weights}{as in \code{\link{glm}}.}
  \item{df.residual}{as in \code{\link{glm}}.}
  \item{df.null}{as in \code{\link{glm}}.}
  \item{y}{as in \code{\link{glm}}.}
  \item{converged}{as in \code{\link{glm}}.}
  \item{boundary}{as in \code{\link{glm}}.}
  \item{ModifiedScores}{the vector of the modified scores for the
    parameters at the final iteration. If \code{pl = TRUE} they are the
    derivatives of the penalized likelihood at the final iteration.}
  \item{FisherInfo}{the Fisher information matrix evaluated at the
    resultant estimates. Only available when \code{method = "brglm.fit"}.}
  \item{hats}{the diagonal elements of the hat matrix. Only available
    when \code{method = "brglm.fit"}}
  \item{nIter}{the number of iterations that were required until
    convergence. Only available when \code{method = "brglm.fit"}.}
  \item{cur.model}{a list with components \code{ar} and \code{at} which
    contains the values of the additive modifications to the responses
    (\code{y}) and to the binomial totals (\code{prior.weights}) at
    the resultant estimates (see \code{\link{modifications}} for more
    information). Only available  when \code{method = "brglm.fit"}.}
  \item{model}{as in \code{\link{glm}}.}
  \item{call}{as in \code{\link{glm}}.}
  \item{formula}{as in \code{\link{glm}}.}
  \item{terms}{as in \code{\link{glm}}.}
  \item{data}{as in \code{\link{glm}}.}
  \item{offset}{as in \code{\link{glm}}.}
  \item{control.glm}{as \code{control} in the result of
    \code{\link{glm}}.}
  \item{control.brglm}{the \code{control.brglm} argument that was passed to
  \code{brglm}. Only available when \code{method = "brglm.fit"}.}
  \item{method}{the method used for fitting the model.}
  \item{contrasts}{as in \code{\link{glm}}.}
  \item{xlevels}{as in \code{\link{glm}}.}
  \item{pl}{logical having the same value with the \code{pl}
    argument passed to \code{brglm}. Only available when \code{method =
      "brglm.fit"}.}
}

\note{
  1. Supported methods for objects of class \code{"brglm"} are:
  \itemize{
    \item{\code{\link{print}}}{through \code{print.brglm}.}
    \item{\code{\link{summary}}}{through \code{summary.brglm}.}
    \item{\code{\link{coefficients}}}{inherited from the
      \code{"glm"} class.}
    \item{\code{\link{vcov}}}{inherited from the\code{"glm"}
      class.}
    \item{\code{\link{predict}}}{inherited from the\code{"glm"}
      class.}
    \item{\code{\link{residuals}}}{inherited from the\code{"glm"}
      class.}
    \item{}{and other methods that apply to objects of class
      \code{"glm"}}
  }
  
  2. A similar implementation of the bias-reduction method could be done
     for every GLM, following Kosmidis (2007) (see also Kosmidis and Firth,
     2009). The full set of families and links will be available in a
     future version. However, bias-reduction is not generally beneficial
     as it is in the binomial family and it could cause inflation of the
     variance (see Firth,  1993).
  
  3. Basically, the differences between maximum likelihood, maximum
     penalized likelihood and the modified scores approach are more
     apparent in small sample sizes, in sparse data sets and in cases
     where complete or quasi-complete separation occurs. Asymptotically
     (as \eqn{n} goes to infinity), the three different approaches are
     equivalent to first order.
         
  4. When an offset is not present in the model, the modified-scores based
     estimates are usually smaller in magnitude than the corresponding
     maximum likelihood estimates, shrinking towards the origin of the
     scale imposed by the link function. Thus, the corresponding
     estimated asymptotic standard errors are also smaller.

     The same is true for the maximum penalized likelihood estimates when
     for example, the logit (where the maximum penalized likelihood and
     modified-scores approaches coincide) or the probit links are
     used. However, generally the maximum penalized likelihood estimates do
     not shrink towards the origin. In terms of mean-value
     parameterization, in the case of maximum penalized likelihood the
     fitted probabilities would shrink towards the point where the Jeffreys
     prior is maximized or equivalently where the  quadratic weights are
     simultaneously maximized (see Kosmidis, 2007).

  5. Implementations of the bias-reduction method for logistic
     regressions can also be found in the \pkg{logistf} package. In
     addition to the obvious advantage of \code{brglm} in the range of
     link functions that can be used (\code{"logit"}, \code{"probit"},
     \code{"cloglog"} and \code{"cauchit"}), \code{brglm} is also more
     efficient computationally. Furthermore, for any user-specified link
     function (see the Example section of \code{\link{family}}), the
     user can specify the corresponding pseudo-data representation to be
     used within \code{brglm} (see \code{\link{modifications}} for
     details).  }

\section{Warnings}{
  1. It is not advised to use methods associated with model comparison
     (\code{\link{add1}}, \code{\link{drop1}},
     \code{\link{anova}}, etc.) on objects of class
     \code{"brglm"}. Model comparison when estimation is performed using
     the modified scores or the penalized likelihood is an on-going
     research topic and will be implemented as soon as it is  concluded.

  2. The use of Akaike's information criterion (AIC) for
     model selection when \code{method = "brglm.fit"} is
     controversial. AIC was developed under the assumptions that (i)
     estimation is by maximum likelihood and (ii) that estimation is
     carried out in a parametric family of distributions that contains
     the \dQuote{true} model. At least the first assumption is not valid
     when using \code{method = "brglm.fit"}. However, since the MLE is
     asymptotically unbiased, asymptotically the modified-scores
     approach is equivalent to maximum likelihood. A more appropriate
     information criterion seems to be Konishi's generalized information
     criterion (see Konishi & Kitagawa, 1996, Sections 3.2 and 3.3), which
     will be implemented in a future version.
  }

\references{
  Bull, S. B., Lewinger, J. B. and Lee, S. S. F. (2007). Confidence
  intervals for multinomial logistic regression in sparse
  data. \emph{Statistics in Medicine} \bold{26}, 903--918.
  
  Firth, D. (1992) Bias reduction, the Jeffreys prior and {GLIM}. In
  \emph{Advances in GLIM and statistical modelling: Proceedings of the
  GLIM 92 conference, Munich}, Eds. L.~Fahrmeir, B.~Francis,
  R.~Gilchrist and G.Tutz, pp. 91--100. New York: Springer.

  Firth, D. (1992) Generalized linear models and Jeffreys priors: An
  iterative generalized least-squares approach. In \emph{Computational
  Statistics I}, Eds. Y. Dodge and J. Whittaker. Heidelberg: Physica-Verlag.

  Firth, D. (1993). Bias reduction of maximum likelihood
  estimates. \emph{Biometrika} \bold{80}, 27--38.
  
  Heinze, G. and Schemper, M. (2002). A solution to the problem of
  separation in logistic regression. \emph{Statistics in Medicine}
  \bold{21}, 2409--2419.
  
  Konishi, S. and Kitagawa, G. (1996). Generalised information criteria
  in model selection. \emph{Biometrika} \bold{83}, 875--890.
  
  Kosmidis, I. (2007). Bias reduction in exponential family nonlinear
  models. \emph{PhD Thesis}, Department of Statistics, University of
  Warwick.

  Kosmidis, I. and Firth, D. (2009). Bias reduction in exponential
  family nonlinear models. \emph{Biometrika}
  \bold{96}, 793--804.
}

\author{Ioannis Kosmidis, \email{ioannis@stats.ucl.ac.uk}}

\seealso{\code{\link{glm}}, \code{\link{glm.fit}}}

\examples{
## Begin Example
data(lizards)
# Fit the GLM using maximum likelihood
lizards.glm <- brglm(cbind(grahami, opalinus) ~ height + diameter +
                  light + time, family = binomial(logit), data=lizards,
                  method = "glm.fit")
# Now the bias-reduced fit:
lizards.brglm <- brglm(cbind(grahami, opalinus) ~ height + diameter +
                  light + time, family = binomial(logit), data=lizards,
                  method = "brglm.fit")
lizards.glm
lizards.brglm
# Other links
update(lizards.brglm, family = binomial(probit))
update(lizards.brglm, family = binomial(cloglog))
update(lizards.brglm, family = binomial(cauchit))
# Using penalized maximum likelihood
update(lizards.brglm, family = binomial(probit), pl = TRUE)
update(lizards.brglm, family = binomial(cloglog), pl = TRUE)
update(lizards.brglm, family = binomial(cauchit), pl = TRUE)
}
\keyword{models}
\keyword{regression}
\keyword{iteration}


