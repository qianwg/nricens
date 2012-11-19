strata.cox.1dim <-
function (time, event, zc, ms.strt) {

  nr   = length(time)
  strt = rep(1, nr)
  nstr = as.integer(nr / ms.strt)

  if (is.null(zc) || nstr == 0 || nstr == 1) return(strt)

  mdl  = coxph(Surv(time, event==1) ~ ., as.data.frame(zc))
  strt = categorize(mdl$linear, nlevel=nstr)

  return(strt)
}
