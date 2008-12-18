\name{summary.etm}
\alias{summary.etm}
\alias{print.summary.etm}
\title{Summary methods for an 'etm' object}
\description{
  Summary method for objects of class \code{etm}
}
\usage{
\S3method{summary}{etm}(object, all = FALSE, ...)
\S3method{print}{summary.etm}(x, ...)
}
\arguments{
  \item{object}{An object of class \code{etm}.}
  \item{all}{If set to TRUE, a data.frame will be computed for all
    transitions that are not 0 in the empirical transition matrix.}
  \item{x}{A \code{summary.cpf} object}
  \item{\dots}{Further arguments}
}
\value{
  A list of data.frames giving the transition probability and stage
  occupation probability estimates. List items are named after the
  possible transition.
  \item{P}{Transtion probability estimates}
  \item{var}{Variance estimates}
  \item{time}{Transition times}
  \item{n.risk}{Number of individuals at risk of experiencing a transition
    just before time t}
  \item{n.event}{Number of events at time t}
}
\author{Arthur Allignol \email{arthur.allignol@fdm.uni-freiburg.de}}
\seealso{\code{\link{etm}}}
\keyword{methods}
\keyword{print}
