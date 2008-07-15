xyplot.etm <- function(x, data=NULL, tr.choice="all", col=1, lty=1,
                       xlab="Time", ylab="Estimated Transition probability", ...) {
  if (!inherits(x, "etm"))
    stop("Argument 'x' must be of class 'etm'")
  ufrom <- unique(x$trans$from)
  uto <- unique(x$trans$to)
  absorb <- setdiff(uto, ufrom)
  nam1 <- dimnames(x$est)[[1]]
  nam2 <- dimnames(x$est)[[2]]
  pos <- c(paste(nam1[!(nam1 %in% as.character(absorb))],
                 nam2[!(nam2 %in% as.character(absorb))]),
           paste(x$trans$from, x$trans$to))
  if (tr.choice[1]=="all") {
    tr.choice <- pos
  }
  else if (sum(tr.choice %in% pos == FALSE) > 0)
    stop("Argument 'tr.choice' and possible transitions must match")
  temp <- sapply(1:length(tr.choice), function(i) {
    strsplit(tr.choice[i], " ")
  })
  temp <- do.call(rbind, temp)
  daten <- lapply(1:length(tr.choice), function(i) {
    est <- x$est[temp[i, 1], temp[i, 2], ]
    cov <- rep(tr.choice[i], length(x$time))
    time <- x$time
    data.frame(est, cov, time)
  })
  daten <- do.call(rbind, daten)
  xyplot(daten$est ~ daten$time | daten$cov, type="s",
         col=col, lty=lty, xlab=xlab, ylab=ylab, ...)
}
  
