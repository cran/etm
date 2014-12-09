\name{ggtransfo.etm}
\alias{ggtransfo.etm}
\alias{ggtransfo}

\title{
  Prepare etm output for plotting with ggplot2
}
\description{
  The \code{ggtransfo} function permits to transform the output of
  \code{etm} such that transition probabilities along with confidence
  intervals can be plotted more easily using the \pkg{ggplot2} package. 
}
\usage{
ggtransfo(x, ...)

\S3method{ggtransfo}{etm}(x, tr.choice, ...)
}

\arguments{
  \item{x}{An object of class 'etm'}
  \item{tr.choice}{Character vector of the form 'c("from to","from
    to")' specifying which transitions should be plotted. Default, all
    the transition probabilities are plotted}
  \item{\dots}{Further arguments. In particular for
  \code{\link{summary.etm}} that is called internally}
}

\value{
  A data frame with the same variables returned by
  \code{\link{summary.etm}}. Addtional variables are
  \item{trans}{transition type. In the same format as given by
  \code{tr.choice}}
  \item{timemax}{Lagged transition times for drawing confidence
  intervals with \code{geom_rect}}

}
\author{
  Arthur Allignol, \email{arthur.allignol@uni-ulm.de}
}

\seealso{
\code{\link{etm}}, \code{\link[ggplot2]{ggplot}},
\code{\link[ggplot2]{geom_rect}} \code{\link[ggplot2]{geom_step}} 
}
\examples{
data(sir.cont)

# Modification for patients entering and leaving a state
# at the same date
# Change on ventilation status is considered
# to happen before end of hospital stay
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
for (i in 2:nrow(sir.cont)) {
  if (sir.cont$id[i]==sir.cont$id[i-1]) {
    if (sir.cont$time[i]==sir.cont$time[i-1]) {
      sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
    }
  }
}

### Computation of the transition probabilities
# Possible transitions.
tra <- matrix(ncol=3,nrow=3,FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE

# etm
tr.prob <- etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)

to_plot <- ggtransfo(tr.prob, tr.choice = c("0 1", "1 0"))

\dontrun{
require(ggplot2)
print(ggplot(to_plot, aes(x = time, y = P)) +
      facet_grid(. ~ trans) +
      geom_step() +
      geom_rect(aes(xmin = time, xmax = timemax, ymin = lower, ymax = upper),
                alpha = 0.5) 
)
}

## abortion
data(abortion)

cif.ab <- etmCIF(Surv(entry, exit, cause != 0) ~ group, abortion,
               etype = cause, failcode = 3)

pcif.ab <- lapply(cif.ab[1:2], ggtransfo,
                  tr.choice = c("0 1", "0 2", "0 3"),
                  ci.fun = "cloglog")

pcif.ab[[1]]$Group <- "Control"
pcif.ab[[2]]$Group <- "Exposed"

pcif.ab <- do.call(rbind, pcif.ab)
pcif.ab$Group <- factor(pcif.ab$Group)
pcif.ab$Group <- relevel(pcif.ab$Group, ref = "Exposed")
pcif.ab$out <- factor(pcif.ab$trans, labels = c("Spontaneous abortion",
                                                 "Induced abortion",
                                                 "Live birth"),
                      ordered = TRUE)

\dontrun{
require(ggplot2)

ggplot(pcif.ab, aes(x = time, y = P)) +
    facet_grid(. ~ out) +
    geom_step(aes(colour = Group), size = 1.3) +
    geom_rect(aes(xmin = time, xmax = timemax, ymin = lower, ymax = upper, fill = Group),
              alpha = 0.5) +
    scale_x_continuous("Week of gestation", limits = c(0, 45)) +
    scale_y_continuous("CIF", limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme(axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          strip.text = element_text(size = 18, face = 2),
          legend.text = element_text(size = 18),
          legend.title = element_blank(),
          legend.position = "top",
          panel.background = element_rect(fill = grey(.93))) +
    scale_colour_brewer(type = "qual", palette = 6) +
    scale_fill_brewer(type = "qual", palette = 6)
}
}

\keyword{hplot}
\keyword{dplot}
