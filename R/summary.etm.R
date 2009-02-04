summary.etm <- function(object, all = FALSE, ...) {
    if (!inherits(object, "etm"))
        stop("'object' must be of class 'etm'")
    if (all) {
        ind <- object$est != 0
        indi <- apply(ind, c(1, 2), function(temp){all(temp == FALSE)})
        tmp <- which(indi == FALSE, arr.ind = TRUE)
        tmp <- tmp[order(tmp[, 1]), ]
        namen <- list(rownames(indi), colnames(indi))
        trs <- lapply(seq_len(NROW(tmp)), function(i) {
            paste(namen[[1]][tmp[i, 1]], namen[[2]][tmp[i, 2]], sep = " ")
        })
        trs <- cbind(trs)
        absorb <- setdiff(levels(object$tran$to), levels(object$trans$from))
        for (i in seq_along(absorb))
            trs <- trs[-grep(paste("^", absorb[i], sep =""), trs, perl = TRUE)]
        trs.sep <- lapply(trs, strsplit, split = " ")
        trs.sep <- matrix(unlist(trs.sep), length(trs.sep), 2, byrow = TRUE)
    }
    else {
        dtrs <- diag(outer(object$state.names, object$state.names, paste))
        absorb <- setdiff(levels(object$tran$to), levels(object$trans$from))
        for (i in seq_along(absorb))
            dtrs <- dtrs[-grep(paste("^", absorb[i], sep =""), dtrs, perl = TRUE)]
        tmp <- paste(object$trans[, 1], object$trans[, 2])
        trs <- c(tmp, dtrs)
        trs.sep <- lapply(trs, strsplit, split = " ")
        trs.sep <- matrix(unlist(trs.sep), length(trs.sep), 2, byrow = TRUE)
    }
    if (!is.null(object$time)) {
        res <- lapply(1:nrow(trs.sep), function(i) {
            P <- object$est[trs.sep[i, 1], trs.sep[i, 2], ]
            time <- object$time
            n.event <- object$n.event[trs.sep[i, 1], trs.sep[i, 2], ]
            if (is.null(dim(object$n.risk)))
                n.risk <- object$n.risk
            else
                n.risk <- object$n.risk[, trs.sep[i, 1]]
            if (is.null(object$cov))
                var <- rep(".", length(P))
            else var <- object$cov[trs[[i]], trs[[i]], ]
            data.frame(P, var, time, n.risk, n.event)
        })
        names(res) <- trs
    }
    else {
        res <- list(P = object$est, s = object$s, t = object$t)
    }
    class(res) <- "summary.etm"
    return(res)
}
