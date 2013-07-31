strata.cox.2dim <-
function (time, event, zc.T, zc.C, m, pca2nd, subs=NULL) {

  nr   = length(time)
  strt = rep(1, nr)
  if (is.null(subs)) subs = rep(TRUE, nr)

  if (is.null(zc.T) && is.null(zc.C)) {
    message("    variables for modeling both T & C ('zc.T','zc.C') are not specified. one stratum.")
    return(list(strt=strt, mdl.T=NULL, mdl.C=NULL, pca=NULL))
  }

  ## subset selection
  time.w  = time[subs]
  event.w = event[subs]
  zc.T.w  = zc.T[subs,]
  zc.C.w  = zc.C[subs,]
  nr.w    = length(time.w)

  nstr = as.integer(nr.w / m)
  if (nstr == 0 || nstr == 1) {
    message("    #of subjects is not enough for building models. one stratum.")
    return(list(strt=strt, mdl.T=NULL, mdl.C=NULL, pca=NULL))
  }

  ## calculating scores for T & C
  mdl.T = try( {coxph(Surv(time.w, event.w==1) ~ ., as.data.frame(zc.T.w))}, silent=TRUE)
  mdl.C = try( {coxph(Surv(time.w, event.w==0) ~ ., as.data.frame(zc.C.w))}, silent=TRUE)

  if (class(mdl.T) == 'try-error' && class(mdl.C) == 'try-error') {
    message("    fail to obtain Cox models for both T & C. one stratum.")
    return(list(strt=strt, mdl.T=NULL, mdl.C=NULL, pca=NULL))
  } else if (class(mdl.T) == 'try-error') {
    if (is.null(zc.T))
      message("    variables for modeling T ('zc.T') is not specified. strata are based on the score for C.")
    else
      message("    fail to obtain Cox model for T. strata are based on the score for C.")
    return(list(strt=strata.1dim(predict(mdl.C, newdata=as.data.frame(zc.C)), m, subs), mdl.T=mdl.T, mdl.C=NULL, pca=NULL))
  } else if (class(mdl.C) == 'try-error') {
    if (is.null(zc.C))
      message("    variables for modeling C ('zc.C') is not specified. strata are based on the score for T.")
    else
      message("    fail to obtain Cox model for C. strata are based on the score for T.")
    return(list(strt=strata.1dim(predict(mdl.T, newdata=as.data.frame(zc.T)), m, subs), mdl.T=NULL, mdl.C=mdl.C, pca=NULL))
  }

  sct = predict(mdl.T, newdata=as.data.frame(zc.T))
  scc = predict(mdl.C, newdata=as.data.frame(zc.C))

  ## fail in model calculations
  if (all(sct == 0) && all(scc == 0)) {
    message("    fail to obtain Cox models for both T & C. one stratum.")
    return(list(strt=strt, mdl.T=NULL, mdl.C=NULL, pca=NULL))
  } else if (all(sct == 0)) {
    message("    fail to obtain Cox model for T. strata are based on the score for T.")
    return(list(strt=strata.1dim(predict(mdl.C, newdata=as.data.frame(zc.C)), m, subs), mdl.T=mdl.T, mdl.C=NULL, pca=NULL))
  } else if (all(scc == 0)) {
    message("    fail to obtain Cox model for C. strata are based on the score for C.")
    return(list(strt=strata.1dim(predict(mdl.T, newdata=as.data.frame(zc.T)), m, subs), mdl.T=NULL, mdl.C=mdl.C, pca=NULL))
  }

  ## pca
  IJ = c(as.integer(nstr / pca2nd), pca2nd)
  if (IJ[1] == 0) IJ[1] = 1

  pca = prcomp(cbind(sct, scc)[subs,], scale=TRUE)

  scr.pca  = cbind(sct/sd(sct), scc/sd(scc)) %*% pca$rotation
  scr.cat1 = categorize(scr.pca[,1], threshold = quantile(scr.pca[subs,1], seq(0, 1, 1/IJ[1])[-1]))
  scr.cat2 = categorize(scr.pca[,2], threshold = quantile(scr.pca[subs,2], seq(0, 1, 1/IJ[2])[-1]))

  pcc = cbind(scr.cat1, scr.cat2)
  unq = unique(pcc)
  unq = unq[order(unq[,1], unq[,2]),]

  for (i in 1:nrow(unq)) {
    ff = apply(pcc, 1, function(xx, uq) { all(xx == uq) }, unq[i,])
    strt[ff] = i
  }
  message('    strata are successfully constructed.')

  return(list(strt=strt, mdl.T=mdl.T, mdl.C=mdl.C, pca=pca))
}
