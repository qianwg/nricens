nribin.count.main <-
function (event, upp, dwn) {

  pr.upp.case = mean(upp & event == 1)
  pr.dwn.case = mean(dwn & event == 1)
  pr.dwn.ctrl = mean(dwn & event == 0)
  pr.upp.ctrl = mean(upp & event == 0)

  nri.case = pr.upp.case - pr.dwn.case
  nri.ctrl = pr.dwn.ctrl - pr.upp.ctrl
  nri      = nri.case + nri.ctrl

  return(c(nri, nri.case, nri.ctrl, pr.upp.case, pr.dwn.case, pr.dwn.ctrl, pr.upp.ctrl))
}
