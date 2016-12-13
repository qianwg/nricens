nricens.km.main <-
function (time, event, upp, dwn, t0) {

  pr.upp = mean(upp)
  pr.dwn = mean(dwn)

  if (pr.upp == 0 && pr.dwn == 0) return(rep(0, 7))

  ## Pr( T>t0 )
  pr.ctrl = get.surv.km(time, as.integer(event==1), t0)

  ## UP  ----  Pr(Up | T<t0), Pr(Up | T>t0)
  pr.ctrl.upp = pr.upp.ctrl = pr.upp.case = 0
  if (pr.upp != 0) {
    pr.ctrl.upp = get.surv.km(time[upp], as.integer(event[upp]==1), t0)
    pr.upp.ctrl = pr.ctrl.upp       * pr.upp / pr.ctrl
    pr.upp.case = (1 - pr.ctrl.upp) * pr.upp / (1 - pr.ctrl)
  }

  ## DOWN  ----  Pr(Down | T<t0), Pr(Down | T>t0)
  pr.ctrl.dwn = pr.dwn.ctrl = pr.dwn.case = 0
  if (pr.dwn != 0) {
    pr.ctrl.dwn = get.surv.km(time[dwn], as.integer(event[dwn]==1), t0)
    pr.dwn.ctrl = pr.ctrl.dwn       * pr.dwn / pr.ctrl
    pr.dwn.case = (1 - pr.ctrl.dwn) * pr.dwn / (1 - pr.ctrl)
  }

  nri.case = pr.upp.case - pr.dwn.case
  nri.ctrl = pr.dwn.ctrl - pr.upp.ctrl
  nri      = nri.case + nri.ctrl

  return(c(nri, nri.case, nri.ctrl, pr.upp.case, pr.dwn.case, pr.dwn.ctrl, pr.upp.ctrl))
}
