#' Computing the stabilized IPCweights
#'
#' @param data a dataframe containing the following variables
#' @param tstart the date of the beginning of the follow-up (in numeric format, with the first being equal at 0)
#' @param tstop the date of the end of the follow-up (in numeric format)
#' @param cens the indicator of treatment censoring (denoted by 1 at the end of the follow-up)
#' @param arm the randomized treatment (2-levels factor)
#' @param bas.cov a vector the baseline covariates
#' @param conf a vector of time-dependent confounders
#' @param trunc an optional fraction for the weights. For instance, when trunc = 0.01,
#' the left tail is truncated to the 1st percentile and the right tail is truncated to the 99th percentile
#'
#' @return the initial dataframe data with stabilized IPCweights as additional arguments. By default, the un-truncated stabilized weights are given. If the trunc option is not NULL then the truncated stabilized weights are also given.
#' @export
#'
#' @importFrom stats quantile
#' @importFrom survival survfit
#' @importFrom survival coxph
#'
#' @references Graffeo, N., Latouche, A., Le Tourneau C., Chevret, S. "An R Package for IPCW: Application to switches in clinical trials" \emph{(submitted)}
#'
#' @examples
#' ## Not run
#' # ipcw(toy.rep, tstart = tstart, tstop = tstop, cens = cens,
#' # arm="arm",
#' # bas.cov = c("age", "arm", "swtrtdt"),
#' # conf = c("TDconf"), trunc = 0.05)
#'
#' # see ?data.ipcw for a complete example
#' @seealso \code{\link{data.ipcw}}
ipcw <- function(data, tstart, tstop, cens, arm, bas.cov, conf, trunc=NULL){

    tempcall <- match.call()

    data$weights <- vector("numeric", length = nrow(data))

    levArms <- levels(data[,arm])

    # ~ baseline covariates in numerator
    num <- paste(c(bas.cov), collapse="+")
    # ~ baseline covariates and time-dependent confoundersin the denominator
    denom <- paste(c(bas.cov, conf), collapse="+")

    # For the first arm
    data1 <- data[data[,arm]==levArms[1],]


    # Cox model for stabilized weights
    #
    fit.cox.num.1 <- coxph(formula = eval(parse(text = paste("Surv(",
                                                             deparse(tempcall$tstart), ", ",
                                                             deparse(tempcall$tstop), ", ",
                                                             deparse(tempcall$cens), ") ~  ",
                                                             num, sep = ""))),
                           data = data1)
    fit.cox.denom.1 <- coxph(eval(parse(text = paste("Surv(",
                                                     deparse(tempcall$tstart), ", ",
                                                     deparse(tempcall$tstop), ", ",
                                                     deparse(tempcall$cens), ") ~  ",
                                                     denom, sep = ""))),
                             data = data1)

    for(i in 1:nrow(data1)) {
        datai <- data1[i,]
        km.num <- survfit(fit.cox.num.1, newdata = datai, type = "kaplan-meier")
        km.denom <- survfit(fit.cox.denom.1, newdata = datai, type = "kaplan-meier")

        data1$weights[i] <- summary(km.num,times=datai$tstart)$surv / summary(km.denom,times=datai$tstart)$surv
    }

    # For the second arm
    data2 <- data[data[,arm]==levArms[2],]
    # Cox model for stabilized weights
    # ~ baseline covariates in numerator
    # ~ baseline covariates and time-dependent confoundersin the denominator
    fit.cox.num.2<- coxph(eval(parse(text = paste("Surv(",
                                                  deparse(tempcall$tstart), ", ",
                                                  deparse(tempcall$tstop), ", ",
                                                  deparse(tempcall$cens), ") ~  ",
                                                  num, sep = ""))),
                          data = data2)
    fit.cox.denom.2 <- coxph(eval(parse(text = paste("Surv(",
                                                     deparse(tempcall$tstart), ", ",
                                                     deparse(tempcall$tstop), ", ",
                                                     deparse(tempcall$cens), ") ~  ",
                                                     denom, sep = ""))),
                               data = data2)

    for(i in 1:nrow(data2)) {
        datai <- data2[i,]
        km.num <- survfit(fit.cox.num.2, newdata = datai, type = "kaplan-meier")
        km.denom <- survfit(fit.cox.denom.2, newdata = datai, type = "kaplan-meier")

        data2$weights[i] <- summary(km.num,times=datai$tstart)$surv / summary(km.denom,times=datai$tstart)$surv
    }


    # rows with 1st level of arm
    rows.1 <- which(data[,arm]==levArms[1])
    data$weights[rows.1] <- data1$weights

    # rows with 2nd level of arm
    rows.2 <- which(data[,arm]==levArms[2])
    data$weights[rows.2] <- data2$weights

    # Truncated weights (optional)
    if (!(is.null(tempcall$trunc))) {
        data$weights.trunc <- data$weights
        data$weights.trunc[data$weights <= quantile(data$weights,
                                                    0 + trunc)] <- quantile(data$weights, 0 +
                                                                                trunc)
        data$weights.trunc[data$weights > quantile(data$weights,
                                                   1 - trunc)] <- quantile(data$weights, 1 -
                                                                               trunc)
    }

    return(data)
}
