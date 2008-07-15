prodint <- function(dna, times, s, t) {
  I <- array(0, dim=dim(dna)[c(1, 2)])
  diag(I) <- 1
  first <- length(times[times <= s]) + 1
  last <- length(times[times <= t])
  if (first >= last) {
    est <- array(I, dim=c(dim(dna)[c(1, 2)], 1))
    time <- NULL
  }
  else {
    est <- array(0, dim=c(dim(dna)[c(1, 2)], (last-first+1)))
    est[, , 1] <- I + dna[, , first]
    j <- 2
    for (i in (first + 1):last) {
      est[, , j] <- est[, , j-1] %*% (I + dna[, , i])
      j <- j + 1
    }
    time <- times[first:last]
  }
  list(est=est, time=time)
}




#####################
### cov of dNA(t) ###
#####################

cov.dNA <- function(nrisk, nev, dd) {
  cov <- matrix(0, dd^2, dd^2)
  for (i in 1:dd) {
    for (j in 1:dd) {
      temp <- matrix(0, dd, dd)
      for (k in 1:dd) {
        for (l in 1:dd) {
          if (k == l) {
            if (k == i) {
              if (l == j) {
                temp[k, l] <- ((nrisk[k] - sum(nev[k, ])) * sum(nev[k, ]))/nrisk[k]^3
              }
              else {
                temp[k, l] <- -(((nrisk[k] - sum(nev[k, ])) * nev[k, j])/nrisk[k]^3)
              }
            }
            else {
              if (i != k & j != k) {
                if (i == j) {
                  temp[k, l] <- ((nrisk[k] - nev[k, i]) * nev[k, i]) / nrisk[k]^3
                }
                else {
                  temp[k, l] <- (-nev[k, i] * nev[k, j]) / nrisk[k]^3
                }
              }
            }
          }
        }
      }
      cov[((i - 1) * dd + 1):(i*dd), ((j - 1) * dd + 1):(j*dd)] <- temp
    }
  }
  cov[is.nan(cov)] <- 0
  for (m in 1:dd^2) {
    for (n in 1:dd^2) {
      if (cov[m, n] != 0) cov[n, m] <- cov[m, n]
    }
  }
  cov
}
      

####################################
### Variance of the AJ estimator ###
####################################

var.aj <- function(est, dna, nrisk, nev, s, t, times) {
  first <- length(times[times <= s]) + 1
  last <- length(times[times <= t])
  if (first >= last) {
    return(0)
  }
  else {
    out <- list()
    cov.dna <- cov.dNA(nrisk[first, ], nev[, , first], dim(nev)[1])
    out[[1]] <- diag(1, dim(nev)[1]^2) %*% cov.dna %*% diag(1, dim(nev)[1]^2)
    d <- dim(nev)[1]
    for (i in 1:length(times[(first + 1):last])) {
      step <- first + i
      cov.dna <- cov.dNA(nrisk[step, ], nev[, , step], d)
      out[[i+1]] <- (t(diag(1, d) + dna[, , step]) %x% diag(1, d)) %*% out[[i]] %*%
        ((diag(1, d) + dna[, , step]) %x% diag(1, d)) +
          (diag(1, d) %x% est[, , i]) %*% cov.dna %*% (diag(1, d) %x% t(est[, , i]))
    }
  }
  return(out)
}

###########
### etm ###
###########      
      
etm <- function(data, state.numbers, tra, cens.name, s, t="last", covariance=TRUE) {
  if (missing(data))
    stop("Argument 'data' is missing with no default")
  if (missing(tra))
    stop("Argument 'tra' is missing with no default")
  if (missing(state.numbers))
    stop("Argument 'state.numbers' is missing with no default")
  if (missing(cens.name))
    stop("Argument 'cens.name' is missing with no default")
  if (missing(s))
    stop("Argument 's' is missing with no default")
  if (!is.data.frame(data))
    stop("Argument 'data' must be a data.frame")
  if (!(xor(sum(c("id", "from", "to", "time") %in% names(data)) != 4,
            sum(c("id", "from", "to", "entry", "exit") %in% names(data)) != 5)))
    stop("'data' must contain the right variables")
  if (nrow(tra) != ncol(tra))
    stop("Argument 'tra' must be a quadratic  matrix.")
  if (sum(diag(tra)) > 0)
    stop("transitions into the same state are not allowed")
  if (nrow(tra) != length(state.numbers)) {
    stop("The row number of 'tra' must be equal to the number of states.")
  }
  if (!is.logical(tra)) {
    stop("'tra' must be a matrix of logical values, which describes the possible transitions.")
  }
  if (length(state.numbers) != length(unique(state.numbers))) {
    stop("The state numbers must be unique.")
  }
  if (!(is.null(cens.name))) {
    if (cens.name %in% state.numbers) {
      stop("The name of the censoring variable just is a name of the model states.")
    }
  }
  ### transitions
  colnames(tra) <- rownames(tra) <- state.numbers
  t.from <- lapply(1:dim(tra)[2], function(i) {
    rep(rownames(tra)[i], sum(tra[i, ]))
  })
  t.from <- unlist(t.from)
  t.to <- lapply(1:dim(tra)[2], function(i) {
    colnames(tra)[tra[i, ]==TRUE]
  })
  t.to <- unlist(t.to)
  trans <- data.frame(from=t.from, to=t.to)
  namen <- paste(trans[, 1], trans[, 2])
  # test on transitions
  test <- unique(paste(data$from, data$to))
  if (!(is.null(cens.name))) {
    ref <- c(paste(trans$from, trans$to), paste(unique(trans$from), cens.name))
  }
  else { ref <- paste(trans$from, trans$to) }
  if (!(all(test %in% ref)==TRUE))
    stop("There is undefined transitions in the data set")
  if (sum(as.character(data$from)==as.character(data$to)) > 0)
    stop("Transitions into the same state are not allowed")
  ### data.frame transformation
  data$from <- as.factor(data$from)
  data$to <- as.factor(data$to)
  if (!(is.null(cens.name))) {
    data$from <- factor(data$from, levels = c(cens.name, state.numbers), ordered = TRUE)
    levels(data$from) <- 0:length(state.numbers)
    data$to <- factor(data$to, levels = c(cens.name, state.numbers), ordered = TRUE)
    levels(data$to) <- 0:length(state.numbers)
  }
  else{
    data$from <- factor(data$from, levels = state.numbers, ordered = TRUE)
    levels(data$from) <- 1:length(state.numbers)
    data$to <- factor(data$to, levels = state.numbers, ordered = TRUE)
    levels(data$to) <- 1:length(state.numbers)
  }
  # if not, put like counting process data
  if ("time" %in% names(data)) {
    data <- data[order(data$id, data$time), ] 
    entree <- double(length(data$time))
    entree[1] <- 0
    for (i in 2:nrow(data)) {
      if (data$id[i] != data$id[i-1]) {
        entree[i] <- 0
      }
      else {
        entree[i] <- data$time[i-1]
      }
    }
    data <- data.frame(id=data$id, from=data$from, to=data$to, entry=entree, exit=data$time)
    if (sum(data$entry < data$exit) != nrow(data))
      stop("Exit time from a state must be > entry time") 
  }
  else {
    if (sum(data$entry < data$exit) != nrow(data))
      stop("Exit time from a state must be > entry time")
  }
  ### Computation of the risk set and dN
  ttime <- c(data$entry, data$exit)
  times <- sort(unique(ttime))
  data$from <- as.integer(as.character(data$from))
  data$to <- as.integer(as.character(data$to))
  temp <- .C("risk_set",
             as.integer(nrow(data)),
             as.integer(length(times)),
             as.integer(c(length(unique(trans[, 1]))), length(times)),
             as.integer(c(dim(tra), length(times))),
             as.double(times),
             as.integer(data$from),
             as.integer(data$to),
             as.double(data$entry),
             as.double(data$exit),
             nrisk=integer(dim(tra)[1] * length(times)),
             ncens=integer(dim(tra)[1] * length(times)),
             nev=integer(dim(tra)[1] * dim(tra)[2] * length(times)),
             PACKAGE="etm")
  nrisk <- matrix(temp$nrisk, ncol=dim(tra)[1], nrow=length(times))
  ncens <- matrix(temp$ncens, ncol=dim(tra)[1], nrow=length(times))
  nev <- array(temp$nev, dim=c(dim(tra), length(times)))
  if (sum(nrisk[1, ]==0)) nrisk[1, ] <- nrisk[2, ]
  ### computation of the increments of the NA estimator
  dna <- array(0, dim=c(dim(tra), length(times)))
  for (tt in 1:length(times)) {
    for (i in 1:dim(nrisk)[2]) {
      if (nrisk[tt, i] != 0) {
        dna[i, , tt] <- nev[i, , tt]/nrisk[tt, i]
        dna[i, i, tt] <- -sum(nev[i, , tt])/nrisk[tt, i]
      }
    }
  }
  ### computation of the Aalen-Johansen estimator
  if (t=="last") t <- times[length(times)]
  if (!(0 <= s & s < t))
    stop("'s' and 't' must be positive, and s < t")
  if (t <= times[1] | s >= times[length(times)])
    stop("'s' or 't' is an invalid time")
  est <- prodint(dna, times, s=s, t=t)
  if (covariance == TRUE) {
    var <- var.aj(est$est, dna, nrisk, nev, s, t, times)
    var <- array(unlist(var), dim=c(dim(nev)[1]^2, dim(nev)[1]^2, length(est$time)))
    pos <- sapply(1:length(state.numbers), function(i) {
      paste(state.numbers, state.numbers[i])
    })
    pos <- matrix(pos)
    dimnames(var) <- list(pos, pos, est$time)
  }
  else var <- NULL
  dimnames(est$est) <- list(state.numbers, state.numbers, est$time)
  dimnames(nev) <- list(state.numbers, state.numbers, times)
  colnames(nrisk) <- state.numbers
  nrisk <- nrisk[, !(colnames(nrisk) %in% setdiff(unique(trans$to), unique(trans$from)))]
  rownames(nrisk) <- times
  res <- list(est=est$est, cov=var, time=est$time, s=s, t=t, trans=trans,
              state.numbers=state.numbers, cens.name=cens.name, nrisk=nrisk, nev=nev)
  class(res) <- "etm"
  res
}
