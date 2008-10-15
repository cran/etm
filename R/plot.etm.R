plot.etm <- function(x, tr.choice, xlab = "Time",
                     ylab = "Transition Probability",
                     legend = TRUE, curvlab, locator = FALSE,
                     coord, col, lty, xlim, ylim, ...) {
    if (!inherits(x, "etm"))
      stop("Argument 'x' must be of class 'etm'")
    if (missing(tr.choice))
        tr.choice <- rownames(x$cov[, , 1])
    if (!all(tr.choice %in% rownames(x$cov[, , 1])))
        stop("Argument 'tr.choice' and the possible transitions must match")
    lt <- length(tr.choice)
    if (missing(col)) {
        col <- rep(1, lt)
    }
    if (missing(lty)) {
        lty <- 1:lt
    }
    if (length(col) < lt)
        col <- col * rep(1, lt)
    if (length(lty) < lt)
        lty <- lty * rep(1, lt)
    if (missing(xlim))
        xlim <- c(x$s, x$t)
    if (missing(ylim))
        ylim <- c(0, 1)
    indices <- do.call(rbind, apply(matrix(tr.choice), 2, strsplit, " ")[[1]])
    probs <- sapply(seq_len(lt), function(i) {
        x$est[indices[i, 1], indices[i, 2], ]
    })
    plot(x$time, probs[, 1], type = "s", col=col[1],
         lty=lty[1], xlim=xlim, ylim=ylim, xlab=xlab,
         ylab=ylab, ...)
    if (lt > 1) {
        for (i in 2:lt) {
            lines(x$time, probs[, i], col=col[i], lty=lty[i],
                  type="s", ...)
        }
    }
    if (legend) {
        if (missing(curvlab)) {
            curvlab <- tr.choice
        }
        args <- list(...)
        i <- pmatch(names(args), names(formals("legend")))
        if (locator) {
            do.call("legend", c(list(locator(1),  curvlab, col=col, lty=lty), args[!is.na(i)]))
        }
        else {
            if (missing(coord)) {
                coord <- c(xlim[1]+0.3,ylim[2]-0.03)
            }
            do.call("legend", c(list(coord[1], coord[2],  curvlab, col=col, lty=lty), args[!is.na(i)]))
        }
    }
}
