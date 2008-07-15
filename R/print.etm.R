print.etm <- function(x, ...) {
  if (!inherits(x, "etm"))
    stop("'x' must be of class 'etm'")
  if (!is.null(x$time)) {
    cat(paste("Time points in the interval (", x$s, ",", x$t, "]", sep=""))
    cat("\n")
    print(x$time)
    cat("\n")
    for (i in 1:length(x$time)) {
      cat(paste("Estimate of P(", x$s, ",", x$time[i], ")", sep=""))
      cat("\n")
      print(x$est[, , i])
      cat("\n")
    }
  }
  else {
    cat(paste("No event in (", x$s, ",", x$t, "]", sep=""))
    cat("\n")
    cat(paste("Estimate of P(", x$s, ",", x$t, "]", sep=""))
    cat("\n")
    print(x$est[, , 1])
    cat("\n")
  }
}
