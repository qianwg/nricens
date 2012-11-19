nricens.km.main <-
function (time, event, upp, dwn, t0, strt, strt.upp, strt.dwn) {

  pr.upp = mean(upp)
  pr.dwn = mean(dwn)

  if (pr.upp == 0 && pr.dwn == 0) return(rep(0, 7))

  ## Pr( T>t0 )
  if (!(length(unique(strt)) == 1 && length(unique(strt.upp)) == 1 && length(unique(strt.dwn)) == 1)) {
    if (pr.upp != 0) strt[upp] = strt.upp[upp] + max(strt)
    if (pr.dwn != 0) strt[dwn] = strt.dwn[dwn] + max(strt)
  }
  pr.ctrl = get.surv.wkm(time, as.integer(event==1), t0, strt)

  ## UP  ----  Pr(Up | T<t0), Pr(Up | T>t0)
  pr.ctrl.upp = pr.upp.ctrl = pr.upp.case = 0
  if (pr.upp != 0) {
    pr.ctrl.upp = get.surv.wkm(time[upp], as.integer(event[upp]==1), t0, strt.upp[upp])
    pr.upp.ctrl = pr.ctrl.upp       * pr.upp / pr.ctrl
    pr.upp.case = (1 - pr.ctrl.upp) * pr.upp / (1 - pr.ctrl)
  }

  ## DOWN  ----  Pr(Down | T<t0), Pr(Down | T>t0)
  pr.ctrl.dwn = pr.dwn.ctrl = pr.dwn.case = 0
  if (pr.dwn != 0) {
    pr.ctrl.dwn = get.surv.wkm(time[dwn], as.integer(event[dwn]==1), t0, strt.dwn[dwn])
    pr.dwn.ctrl = pr.ctrl.dwn       * pr.dwn / pr.ctrl
    pr.dwn.case = (1 - pr.ctrl.dwn) * pr.dwn / (1 - pr.ctrl)
  }

  nri.case = pr.upp.case - pr.dwn.case
  nri.ctrl = pr.dwn.ctrl - pr.upp.ctrl
  nri      = nri.case + nri.ctrl

  return(c(nri, nri.case, nri.ctrl, pr.upp.case, pr.dwn.case, pr.dwn.ctrl, pr.upp.ctrl))
}
