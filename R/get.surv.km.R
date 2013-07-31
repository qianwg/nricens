get.surv.km <-
function(time, event, t0, subs=NULL) {

  if (!is.null(subs)) {
    if (sum(subs) == 0)
      stop("error in get.surv.km")

    time  = time[subs]
    event = event[subs]
  }

  svf  = survfit(Surv(time, event) ~ 1, se.fit=FALSE)
  surv = stepfun(svf$time, c(1, svf$surv))(t0)

  return(surv)
}
