\name{print-summary}
\alias{print.etm}
\alias{summary.etm}
\title{Print and summary methods for object of class 'etm'}
\description{
  Print and summary methods for objects of class \code{etm}.
}
\usage{
\S3method{print}{etm}(x, ...)
\S3method{summary}{etm}(object, ...)
}
\arguments{
  \item{x}{An object of class \code{etm}.}
  \item{object}{An object of class \code{etm}.}
  \item{\dots}{Further arguments for print or summary.}
}
\details{
  The print method displays the estimated matrix of transition
  probabilities for each event times in \eqn{(s, t]}{(s, t]}.
  
  The summary methods provides a description of the multistate model,
  along with the estimation of \eqn{P(s, t]}{P(s, t]} and cov\eqn{(P(s,
  t])}{(P(s, t])}.
}
\value{
  No value returned
}
\author{Arthur Allignol, \email{arthur.allignol@fdm.uni-freiburg.de}}
\seealso{\code{\link{etm}}}
\keyword{print}
