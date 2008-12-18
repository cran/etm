plot.etm <- function(x, tr.choice = "all", xlab = "Time",
                     ylab = "Transition Probability",
                     col, lty, xlim, ylim, legend = TRUE,
                     legend.pos, curvlab, legend.bty = "n", ...) {
    if (!inherits(x, "etm"))
        stop("'x' must be a 'etm' object")
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
    else if (sum(tr.choice %in% rownames(x$cov) == FALSE) > 0)
        stop("Argument 'tr.choice' and possible transitions must match")
    lt <- length(tr.choice)
    indices <- do.call(rbind, apply(matrix(tr.choice), 2, strsplit, " ")[[1]])
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
    probs <- sapply(seq_len(lt), function(i) {
        x$est[indices[i, 1], indices[i, 2], ]
    })
    plot(x$time, x$time, xlim = xlim, ylim = ylim,
         type = "n", xlab = xlab, ylab = ylab, ...)
    for (i in seq_along(tr.choice)) {
        lines(x$time, probs[, i], col = col[i], lty = lty[i], type = "s", ...)
    }
    if (legend) {
        if (missing(legend.pos))
            legend.pos <- "topright"
        if (missing(curvlab))
            curvlab <- tr.choice
        if (is.list(legend.pos)) legend.pos <- unlist(legend.pos)
        if (length(legend.pos) == 1) {
            xx <- legend.pos
            yy <- NULL
        }
        if (length(legend.pos) == 2) {
            xx <- legend.pos[1]
            yy <- legend.pos[2]
        }
        args <- list(...)
        ii <- pmatch(names(args),
                     names(formals("legend")[-charmatch("bty",names(formals("legend")))]))
        do.call("legend", c(list(xx, yy, curvlab, col=col, lty=lty, bty = legend.bty),
                            args[!is.na(ii)]))
    }
    invisible(list(s = x$s, t = x$t))
}
