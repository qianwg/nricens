strata.cox.2dim <-
function (time, event, zcT, zcC, ms.strt, pca2nd, subs=NULL) {

  nr   = length(time)
  strt = rep(1, nr)
  if (is.null(zcT)) return(strt)

  if (is.null(subs)) subs = rep(TRUE, nr)
  if (is.null(zcC))  zcC  = zcT

  ## samples for scores in PCA
  time.w  = time[subs]
  event.w = event[subs]
  zcT.w   = zcT[subs,]
  zcC.w   = zcC[subs,]
  nr.w    = length(time.w)

  nstr = as.integer(nr.w / ms.strt)
  if (nstr == 0 || nstr == 1) return(strt)

  ## calculating scores
  mdlT = try( {coxph(Surv(time.w, event.w==1) ~ ., as.data.frame(zcT.w))}, silent=T)
  mdlC = try( {coxph(Surv(time.w, event.w==0) ~ ., as.data.frame(zcC.w))}, silent=T)

  if (class(mdlT) == 'try-error' && class(mdlT) == 'try-error') {
    return(strt)
  } else if (class(mdlT) == 'try-error') {
    return(strata.cox.1dim(time, as.integer(event==0), zcC, ms.strt))
  } else if (class(mdlC) == 'try-error') {
    return(strata.cox.1dim(time, as.integer(event==1), zcT, ms.strt))
  }

  sct = predict(mdlT, newdata=as.data.frame(zcT))
  scc = predict(mdlC, newdata=as.data.frame(zcC))

  ## fail in model calculations
  if (all(sct == 0) && all(scc == 0)) {
    return(strt)
  } else if (all(sct == 0)) {
    thr  = quantile(scc[subs], seq(0, 1, 1/nstr)[-1])
    strt = categorize(scc, threshold=thr)
    return(strt)
  } else if (all(scc == 0)) {
    thr  = quantile(sct[subs], seq(0, 1, 1/nstr)[-1])
    strt = categorize(sct, threshold=thr)
    return(strt)
  }

  ## pca
  IJ = c(as.integer(nstr / pca2nd), pca2nd)
  if (IJ[1] == 0) IJ[1] = 1

  pca = prcomp(cbind(sct, scc)[subs,], scale=TRUE)

  scr.pca = cbind(sct/sd(sct), scc/sd(scc)) %*% pca$rotation
  thr1    = quantile(scr.pca[subs,1], seq(0, 1, 1/IJ[1])[-1])
  thr2    = quantile(scr.pca[subs,2], seq(0, 1, 1/IJ[2])[-1])

  scr.cat1 = categorize(scr.pca[,1], threshold=thr1)
  scr.cat2 = categorize(scr.pca[,2], threshold=thr2)

  pcc = cbind(scr.cat1, scr.cat2)
  unq = unique(pcc)
  unq = unq[order(unq[,1], unq[,2]),]

  for (i in 1:nrow(unq)) {
    ff = apply(pcc, 1, function(xx, uq) { all(xx == uq) }, unq[i,])
    strt[ff] = i
  }

  return(strt)
}
