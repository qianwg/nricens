get.surv.km <-
function(time, event, t0, subst=NULL) {

  if (!is.null(subst)) {
    if (sum(subst) == 0) {
      print("error in get.surv.km")
      return(1)
    }

    time  = time[subst]
    event = event[subst]
  }

  svf  = survfit(Surv(time, event) ~ 1, se.fit=F)
  surv = stepfun(svf$time, c(1, svf$surv))(t0)

  return(surv)
}
