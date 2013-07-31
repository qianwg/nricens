get.surv.wkm <-
function (time, event, t0, strt) {

  surv = 0
  for (i in unique(strt)) {
    ff = strt == i
    if (mean(ff) != 0)
      surv = surv + get.surv.km(time, event, t0, ff) * mean(ff)
  }

  return(surv)
}
