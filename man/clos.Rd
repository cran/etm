\name{clos}
\alias{clos}
\title{Change in Length of Stay}
\description{
  The function estimates the expected change in length of stay (LOS)
  associated with an intermediate event in the same way as does the
  \code{\link[changeLOS]{clos}} function in the \pkg{changeLOS}
  package. The difference between the 2 implementations is that this
  function can handle left-truncated data, and does not require a
  competing outcome.
}
\usage{
clos(x, aw = FALSE)
}
\arguments{
  \item{x}{An object of class \code{etm}. Argument \code{delta.na} in
    \code{\link{etm}} must be set to \code{TRUE} in order to use this
    function.}
  \item{aw}{Logical. Whether to compute the expected change of LOS using
    alternative weighting. Default is \code{FALSE}.}
}
\details{
  See \code{\link[changeLOS]{clos}} for more details on the method
}
\value{
  An object of class \code{clos.etm} with the following components:
  \item{e.phi}{Change in length of stay}
  \item{phi.case}{Estimates of \eqn{E(\mbox{LOS} | X_s =
      \mbox{intermediate event})}{E(LOS | X_s = intermediate event)} for
    all observed transition times \eqn{s}{s}, where
    \eqn{X_s}{X_s}denotes the state by time \eqn{s}{s}}
  \item{phi.control}{Estimates of \eqn{E(\mbox{LOS} | X_s =
      \mbox{initial state})}{E(LOS|X_s = initial state)} for
    all observed transition times \eqn{s}{s}.}
  \item{e.phi2}{Weighted average of the difference between
    \code{phi2.case} and \code{phi2.control}.}
  \item{phi2.case}{Estimates of \eqn{E(\mbox{LOS}
      \mathbf{1}(X_{\mbox{LOS}} = \mbox{discharge}) | X_s =
      \mbox{intermediate event})}{E(LOS \strong{1}(X_LOS = discharge)|X_s =
      intermediate event)}, where \eqn{\mathbf{1}}{\strong{1}} denotes
    the indicator function.}
  \item{phi2.control}{\eqn{E(\mbox{LOS}
      \mathbf{1}(X_{\mbox{LOS}} = \mbox{discharge}) | X_s =
      \mbox{initial state})}{E(LOS \strong{1}(X_LOS = discharge)|X_s =
      initial state)}.}
  \item{e.phi3}{Weighted average of the difference between
    \code{phi3.case} and \code{phi3.control}.}
  \item{phi3.case}{Estimates of \eqn{E(\mbox{LOS}
      \mathbf{1}(X_{\mbox{LOS}} = \mbox{death}) | X_s =
      \mbox{intermediate event})}{E(LOS \strong{1}(X_LOS = death)|X_s =
      intermediate event)}.}
  \item{phi3.control}{\eqn{E(\mbox{LOS}
      \mathbf{1}(X_{\mbox{LOS}} = \mbox{death}) | X_s =
      \mbox{initial state})}{E(LOS \strong{1}(X_LOS = death)|X_s =
      initial state)}.}
  \item{weights}{Weights used to compute the weighted averages.}
  \item{w.time}{Times at which the weights are computed.}
  \item{time}{All transition times.}
  \item{e.phi.weights.1}{Expected change in LOS using \code{weights.1}}
  \item{e.phi.weights.other}{Expected change in LOS using
    \code{weights.other}}
  \item{weights.1}{Weights corresponding to the conditional waiting
    time in the intial state given one experiences the intermediate event.}
  \item{weights.other}{Weights corresponding to the conditional
    waiting time given one does not experience the intermediate event.} 
}
\author{Arthur Allignol \email{arthur.allignol@fdm.uni-freiburg.de},
  Matthias Wangler, Jan Beyersmann}
\seealso{\code{\link[changeLOS]{clos}}, \code{\link{etm}}}
\examples{
require(changeLOS)
data(los.data) # in package changeLOS

## putting los.data in the long format (see changeLOS)
my.observ <- prepare.los.data(x=los.data)

tra <- matrix(FALSE, 4, 4)
tra[1, 2:4] <- TRUE
tra[2, 3:4] <- TRUE

tr.prob <- etm(my.observ, c("0","1","2","3"), tra, NULL, 0)

cLOS <- etm::clos(tr.prob)
plot(cLOS)


### Compute bootstrapped SE

## function that performs the bootstrap
## nboot: number of bootstrap samples. Other arguments are as in etm()
boot.clos <- function(data, state.names, tra, cens.name, s = 0, nboot) {
    res <- double(nboot)
    for (i in seq_len(nboot)) {
        index <- sample(unique(data$id), replace = TRUE)
        inds <- new.id <- NULL
        for (j in seq_along(index)){
            ind <- which(data$id == index[j])
            new.id <- c(new.id, rep(j, length(ind)))
            inds <- c(inds, ind)
        }
        dboot <- cbind(data[inds, ], new.id)
        dboot[, which(names(dboot) == "id")]
        dboot$id <- dboot$new.id
        tr.prob <- etm(dboot, state.names, tra, cens.name, s, cova = FALSE)
        res[i] <- etm::clos(tr.prob)$e.phi
    }
    res
}

## bootstrap
se <- sqrt(var(boot.clos(my.observ, c("0","1","2","3"), tra, NULL, 0,
                         nboot = 10)))
}

\keyword{survival}
