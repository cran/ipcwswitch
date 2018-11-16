#' Function to replicate the rows so that each patients' follow-up is split
#' according to all event times (times parameter)
#' up to each patient's end time
#'
#' @param data a dataframe containing the following variables
#' @param tstart the date of the beginning of the follow-up (in numeric format, with the first being equal at 0)
#' @param tstop the date of the end of the follow-up (in numeric format)
#' @param event the indicator of failure (a death is denoted by 1 at the end of the follow-up)
#' @param cens the indicator of treatment censoring (denoted by 1 at the end of the follow-up)
#' @param times a vector of times (in numeric format) inidicating the times according to which
#' the rows have to be split
#'
#' @return a formatted dataframe with the rows replicated according to the provided times parameter
#' @export
#'
#' @importFrom survival survSplit
#'
#' @references Graffeo, N., Latouche, A., Le Tourneau C., Chevret, S. "An R Package for IPCW: Application to switches in clinical trials" \emph{(submitted)}
#'
#'
#' @examples
#' # To obtain the times parameter, we can apply the timesTokeep function on the same
#' # dataframe in the wide format
#' kept.t <- timesTokeep(toydata, id = "id",
#' tstart = "randt", tstop = "lastdt",
#' mes.cov = list(c("ps1", "ps2", "ps3")),
#' time.cov = list(c("randt", "dt2", "dt3")))
#' # Now, we can build the long format
#' toy.long <- wideToLongTDC(data = toydata, id = "id",
#' tstart = "randt", tstop = "lastdt", event = "status",
#' bas.cov = c("age", "arm", "swtrtdt"),
#' mes.cov = list(TDconf = c("ps1", "ps2", "ps3")),
#' time.cov = list(c("randt", "dt2", "dt3")),
#' times = kept.t[[1]])
#' # Put dates in numeric format with tstart at 0
#' toy.long$tstart <- as.numeric(toy.long$tstart)
#' toy.long$tstop <- as.numeric(toy.long$tstop)
#' toy.long$swtrtdt <- as.numeric(toy.long$swtrtdt)
#' tabi <- split(toy.long, toy.long$id)
#' L.tabi   <- length(tabi)
#' tablist <- lapply(1:L.tabi, function(i){
#'     refstart <- tabi[[i]]$tstart[1]
#'     tabi[[i]]$tstart  <- tabi[[i]]$tstart - refstart
#'     tabi[[i]]$tstop <- tabi[[i]]$tstop - refstart
#'     tabi[[i]]$swtrtdt <- tabi[[i]]$swtrtdt - refstart
#'     return(tabi[[i]])
#'     })
#'     toy.long <- do.call( rbind, tablist )
#' # Patients are censored when initiating the other arm treatment, that is, at time swtrtdt
#' toy.long2 <- cens.ipw(toy.long, id = "id", tstart = "tstart", tstop = "tstop",
#' event = "event", arm = "arm",
#' realtrt = FALSE, censTime ="swtrtdt")
#' # We collect all event times (death and treatment censoring)
#' rep.times <- unique(c(toy.long2$tstop[toy.long2$cens==1],
#' toy.long2$tstop[toy.long2$event==1]))
#' # Now, we can replicate the rows
#' toy.rep   <- replicRows(toy.long2, tstart = "tstart", tstop = "tstop",
#'                         event = "event", cens = "cens", times = rep.times)
#' toy.rep
#' @seealso \code{\link{cens.ipw}}, \code{\link{data.ipcw}}, \code{\link{timesTokeep}}, \code{\link{wideToLongTDC}}
replicRows <- function(data, tstart, tstop, event, cens, times){


  # gathering all unique times
  alltimes <- sort(unique(unlist(times)))


  tabi <- split(data, data[,"id"])
  L.tabi   <- length(tabi)
  tablist <- lapply(1:L.tabi, function(i){
    l.tabi <- nrow(tabi[[i]])
    d1 <- list(); d2 <- list()
    for(j in seq(l.tabi)){
      cut.t <- alltimes[(alltimes > tabi[[i]][j,tstart]) & (alltimes < tabi[[i]][j,tstop])]
      d1[[j]] <- survSplit(tabi[[i]][j,],
                      cut = cut.t,
                      end = tstop,
                      event = event)
      d2[[j]] <- survSplit(tabi[[i]][j,],
                      cut = cut.t,
                      end = tstop,
                      event = cens)
      d1[[j]]$cens <- d2[[j]]$cens
    }

    return(do.call(rbind, d1))
  })
  SHIlong2 <- do.call( rbind, tablist )

}
