nricens.ipw.main <-
function (time, event, upp, dwn, t0, strt) {

  case = time <= t0 & event == 1
  ctrl = time >  t0

  w.case = rep(NA, length(time))
  w.ctrl = rep(NA, length(time))

  unq = unique(strt)
  for (i in unq) {
    ff = strt == i
    wk = get.surv.km(time[ff], as.integer(event[ff]==0), c(time[ff], t0))
    w.case[ff] = 1 / wk[-length(wk)]
    w.ctrl[ff] = 1 / wk[length(wk)]
  }

  w.case[!case] = 0
  w.ctrl[!ctrl] = 0
  
  pr.upp.case = sum(w.case * upp) / sum(w.case)
  pr.dwn.case = sum(w.case * dwn) / sum(w.case)
  pr.dwn.ctrl = sum(w.ctrl * dwn) / sum(w.ctrl)
  pr.upp.ctrl = sum(w.ctrl * upp) / sum(w.ctrl)

  nri.case = pr.upp.case - pr.dwn.case
  nri.ctrl = pr.dwn.ctrl - pr.upp.ctrl
  nri      = nri.case + nri.ctrl

  return(c(nri, nri.case, nri.ctrl, pr.upp.case, pr.dwn.case, pr.dwn.ctrl, pr.upp.ctrl))
}
