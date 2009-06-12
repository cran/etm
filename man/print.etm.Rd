\name{print.etm}
\alias{print.etm}
\title{Print method for object of class 'etm'}
\description{
  Print method for objects of class \code{etm}.
}
\usage{
\S3method{print}{etm}(x, covariance = TRUE, whole = TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{etm}.}
  \item{covariance}{Whether print the covariance matrix. Default is
    TRUE}
  \item{whole}{Whether to plot the entire covariance matrix. If set to
    FALSE, rows and columns containing only 0 will be removed for
    printing.}
  \item{\dots}{Further arguments for print or summary.}
}
\details{
  The function prints a matrix giving the possible transitions, along
  with the estimates of \eqn{P(s, t)}{P(s, t)} and \eqn{cov(P(s,
    t))}{cov(P(s, t))}.
}
\value{
  \code{x}, with the invisible flag set.
}
\author{Arthur Allignol, \email{arthur.allignol@fdm.uni-freiburg.de}}
\seealso{\code{\link{etm}}}
\keyword{print}