nricens <-
function (time, event, z.std, z.new, t0, cut, niter=1000, zcT=NULL, zcC=NULL, ms.strt=30, pca2nd=1,
          region.method='cutoff', point.method='km', alpha=0.05) {

  est = rep(NA, 7)
  ci  = rep(NA, 14)

  ##
  ## DH & DL
  wk = get.uppdwn(time, event, z.std, z.new, t0, cut, region.method, msg=FALSE)
  upp = wk[[1]]
  dwn = wk[[2]]

  ##
  ## point estimation
  if (point.method=='ipw') {
    strt = strata.cox.1dim(time, as.integer(event==0), zcC, ms.strt)
    est  = nricens.ipw.main(time, event, upp, dwn, t0, strt)
    strt.upp = strt.dwn = strt

  } else if (point.method=='km') {
    strt.upp = strata.cox.2dim(time, event, zcT, zcC, ms.strt, pca2nd, upp)
    strt.dwn = strata.cox.2dim(time, event, zcT, zcC, ms.strt, pca2nd, dwn)
    strt     = strata.cox.2dim(time, event, zcT, zcC, ms.strt, pca2nd, !upp & !dwn)
    est = nricens.km.main(time, event, upp, dwn, t0, strt, strt.upp, strt.dwn)
  }

  ##
  ## interval estimation
  if (niter > 0) {
    N    = length(time)
    samp = matrix(NA, niter, 7)

    if (point.method == 'km') {
      for (b in 1:niter) {
        f   = as.integer(runif(N, 0, N)) + 1
        wk  = get.uppdwn(time[f], event[f], z.std[f,], z.new[f,], t0, cut, region.method, msg=FALSE)
        upp = wk[[1]]
        dwn = wk[[2]]
        samp[b,] = nricens.km.main(time[f], event[f], upp, dwn, t0, strt[f], strt.upp[f], strt.dwn[f])
      }

    } else if (point.method == 'ipw') {
      for (b in 1:niter) {
        f   = as.integer(runif(N, 0, N)) + 1
        wk  = get.uppdwn(time[f], event[f], z.std[f,], z.new[f,], t0, cut, region.method, msg=FALSE)
        upp = wk[[1]]
        dwn = wk[[2]]
        samp[b,] = nricens.ipw.main(time[f], event[f], upp, dwn, t0, strt[f])
      }
    }

    ci = as.numeric(apply(samp, 2, quantile, c(alpha/2, 1-alpha/2), na.rm=T, type=2))
  }


  result = as.data.frame(cbind(est, matrix(ci, ncol=2, byrow=T)))
  names(result) = c('Estimate', 'Lower', 'Upper')
  row.names(result) = c('NRI','NRI+','NRI-','Pr(Up|Case)','Pr(Down|Case)','Pr(Down|Ctrl)','Pr(Up|Ctrl)')

  return(result)
}
