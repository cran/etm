\name{xyplot.etm}
\alias{xyplot.etm}
\title{xyplot method for object of class 'etm'}
\description{
  xyplot function for objects of class \code{etm}. Estimates of the 
  transition probabilities are plotted as a function of time for all the
  transitions specified by the user.
}
\usage{
\S3method{xyplot}{etm}(x, data = NULL, tr.choice = "all", col = 1, lty = 1, xlab =
"Time", ylab = "Estimated Transition probability", ...)
}
\arguments{
  \item{x}{An object of class \code{etm}.}
  \item{data}{\emph{Useless}.}
  \item{tr.choice}{A character vector of the form c("from to", "from
    to", ...) specifying the transition probabilities to be plotted. By
    default, all the transition probabilities are displayed.}
  \item{col}{Color of the curve.}
  \item{lty}{Type of the line.}
  \item{xlab}{x-axis label. Default is "Time".}
  \item{ylab}{y-axis label. Default is "Estimated transition probability".}
  \item{\dots}{Further arguments for \code{xyplot}.}
}
\value{
  An object of class \code{trellis}.
}
\author{Arthur Allignol, \email{arthur.allignol@fdm.uni-freiburg.de}}
\seealso{\code{\link{etm}}, \code{\link[lattice]{xyplot}}}
\keyword{device}
