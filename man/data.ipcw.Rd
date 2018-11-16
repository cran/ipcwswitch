\name{data.ipcw}
\alias{data.ipcw}
\docType{data}
\title{
A not so short example dataset
}
\description{
Dataset \code{data.ipcw} mimics randomized clinical trials data with 100 patients. It has two parallel arms with a repeated measurement of a time-varying binary covariate, which could be the time-varying confounder acting both on the survival and treatment censoring. The main difference wrt \code{toydata} is that there is enough information to run a complete analysis.
}
\usage{data("data.ipcw")}
\format{
  A data frame with 100 observations on the following 14 variables.
  \describe{
    \item{\code{id}}{a numeric vector corresponding to the patient's identifier    }
    \item{\code{randt}}{a Date corresponding to the date of the randomization     visit}
    \item{\code{lastdt}}{a Date corresponding to the date of latest news}
    \item{\code{status}}{a numeric vector The value equals to 1 if the patient     dies at lastdt (and 0 otherwise)}
    \item{\code{age}}{a numeric vector containing patient’s age (in years) at     randomization}
    \item{\code{score1}}{a numeric vector of the repeated measurement
        happening on date \code{randt}. This first value equals 0.}
    \item{\code{score2}}{a numeric vector of the repeated measurement
        happening on date \code{date.mes2}. This second value equals 1.}
    \item{\code{score3}}{a numeric vector of the repeated measurement
        happening on date \code{date.mes3}. This second value equals 0.}
    \item{\code{score4}}{a numeric vector of the repeated measurement
        happening on date \code{date.mes4}. This second value equals 0.}
    \item{\code{date.mes2}}{the Date of measurement of \code{score2}.}
    \item{\code{date.mes3}}{the Date of measurement of \code{score3}.}
    \item{\code{date.mes4}}{the Date of measurement of \code{score4}.}
    \item{\code{arm}}{a factor with levels \code{trt} \code{pcb}}
    \item{\code{swtrtdt}}{the Date when the patient initiates
        the other arm treatment}
  }
}
\examples{
\donttest{
# To obtain the times parameter, we can apply the timesTokeep function on the same
# dataframe in the wide format
kept.t <- timesTokeep(data.ipcw, id = "id",
                      tstart = "randt", tstop = "lastdt",
                      mes.cov = list(c("score1", "score2", "score3", "score4")),
                      time.cov = list(c("randt", "date.mes2", "date.mes3", "date.mes4")))
# Now, we can build the long format
toy.long <- wideToLongTDC(data = data.ipcw, id = "id",
                          tstart = "randt", tstop = "lastdt", event = "status",
                          bas.cov = c("age", "arm", "swtrtdt"),
                          mes.cov = list(score = c("score1", "score2", "score3", "score4")),
                          time.cov = list(c("randt", "date.mes2", "date.mes3", "date.mes4")),
                          times = kept.t[[1]])
# Put dates in numeric format with tstart at 0
toy.long$tstart <- as.numeric(toy.long$tstart)
toy.long$tstop <- as.numeric(toy.long$tstop)
toy.long$swtrtdt <- as.numeric(toy.long$swtrtdt)

tabi <- split(toy.long, toy.long$id)
L.tabi   <- length(tabi)
tablist <- lapply(1:L.tabi, function(i){
    refstart <- tabi[[i]]$tstart[1]
    tabi[[i]]$tstart  <- tabi[[i]]$tstart - refstart
    tabi[[i]]$tstop <- tabi[[i]]$tstop - refstart
    tabi[[i]]$swtrtdt <- tabi[[i]]$swtrtdt - refstart
    return(tabi[[i]])
})

toy.long <- do.call( rbind, tablist )

# Patients are censored when initiating the other arm treatment, that is, at time swtrtdt
toy.long2 <- cens.ipw(toy.long, id = "id", tstart = "tstart", tstop = "tstop",
                      event = "event", arm = "arm",
                      realtrt = FALSE, censTime ="swtrtdt")
# We collect all event times (death and treatment censoring)
rep.times <- unique(c(toy.long2$tstop[toy.long2$cens==1],
                      toy.long2$tstop[toy.long2$event==1]))
# Now, we can replicate the rows
toy.rep   <- replicRows(toy.long2, tstart = "tstart", tstop = "tstop",
                        event = "event", cens = "cens", times = rep.times)

# Estimation of the stabilized weights
library(survival)
toy.res <- ipcw(toy.rep, tstart = tstart, tstop = tstop, cens = cens,
                arm="arm",
                bas.cov = c("age", "arm"),
                conf = c("score"), trunc = 0.05)

# Using the IPCW weights in Cox likelihood...
fit.stab.w <- coxph(Surv(tstart, tstop, event) ~ arm + age
                    + cluster(id),
                    data = toy.res, weights = toy.res$weights.trunc)
fit.stab.w

}
}
\seealso{
\code{\link{cens.ipw}}, \code{\link{ipcw}}, \code{\link{replicRows}}, \code{\link{timesTokeep}}, \code{\link{wideToLongTDC}}
}
\keyword{datasets}
