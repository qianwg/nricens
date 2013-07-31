nricens <-
function (time=NULL, event=NULL, mdl.std=NULL, mdl.new=NULL, z.std=NULL, z.new=NULL, p.std=NULL, p.new=NULL,
  t0=NULL, updown='category', cut=NULL, point.method='km', zc.T=NULL, zc.C=NULL, mdl.C=NULL, m=100, pca2nd=1,
  niter=1000, alpha=0.05, msg=TRUE) {

  ##
  ## type of calculation
  flag.mdl = !is.null(mdl.std) && !is.null(mdl.new)
  flag.prd = !is.null(z.std)   && !is.null(z.new)
  flag.rsk = !is.null(p.std)   && !is.null(p.new)

  ##
  ## check standard & new model
  if (flag.mdl) {
    if (is.null(time))   time  = as.numeric(mdl.std$y[,1])
    if (is.null(event))  event = as.numeric(mdl.std$y[,2])

    if (class(mdl.std)=='survreg')
      time = survreg.distributions[[mdl.std$dist]]$itrans(time)

    if (is.null(mdl.std$x) || is.null(mdl.new$x))
      stop("\n\nmodel object does not contain predictors. pls set x=TRUE for model calculation.\n\n")

    z.std = mdl.std$x
    z.new = mdl.new$x
    if (any(class(mdl.std) == 'coxph')) {
      mdl.std = coxph(Surv(time, event) ~ ., data=as.data.frame(cbind(time, event, z.std)))
      mdl.new = coxph(Surv(time, event) ~ ., data=as.data.frame(cbind(time, event, z.new)))
    } else if (any(class(mdl.std) == 'survreg')) {
      mdl.std = survreg(Surv(time, event) ~ ., dist=mdl.std$dist, data=as.data.frame(cbind(time, event, z.std)))
      mdl.new = survreg(Surv(time, event) ~ ., dist=mdl.new$dist, data=as.data.frame(cbind(time, event, z.new)))
    }

  } else if (flag.prd) {
    mdl.std = coxph(Surv(time, event==1) ~ ., as.data.frame(cbind(time, event, z.std)))
    mdl.new = coxph(Surv(time, event==1) ~ ., as.data.frame(cbind(time, event, z.new)))
    message("\nSTANDARD prediction model (Cox model):")
    print(summary(mdl.std)$coef)
    message("\nNEW prediction model (Cox model):")
    print(summary(mdl.new)$coef)

  } else if (!flag.mdl && !flag.prd && !flag.rsk) {
    stop("\n\neither one of 'time, event, z.std, z.new', 'time, event, p.std, p.new', and 'mdl.std, mdl.new' should be specified.\n\n")
  }

  if (is.null(cut))
    stop("\n\n'cut' is empty")

  get.risk = NULL
  if (any(class(mdl.std) == 'coxph'))
    get.risk = get.risk.coxph
  else if (any(class(mdl.std) == 'survreg'))
    get.risk = get.risk.survreg
  else if (!flag.rsk)
    stop("'coxph' and 'survreg' are allowed for risk prediction models.")

  objs = list(mdl.std, mdl.new, z.std, z.new, p.std, p.new)

  ##
  ## DH & DL
  wk  = get.uppdwn(time, event, objs, flag.mdl, flag.prd, flag.rsk, t0, updown, cut, get.risk, msg=msg)
  upp = wk[[1]]
  dwn = wk[[2]]
  ret = list(mdl.std=mdl.std, mdl.new=mdl.new, p.std=wk[[3]], p.new=wk[[4]], up=upp, down=dwn, rtab=wk[[5]], rtab.case=wk[[6]], rtab.ctrl=wk[[7]])

  ##
  ## point estimation
  est = rep(NA, 7)
  ci  = rep(NA, 14)
  if (point.method=='ipw' || point.method=='modipw') {

    if (point.method=='modipw' && is.null(mdl.C) && is.null(zc.C)) {
      message("\nNRI estimation by modified IPW estimator:")
      message("  neither model nor predictor for censoring time is specified. standard IPW estimator is used.")
      strt = rep(1, length(time))

    } else if (point.method=='modipw') {
      message("\nNRI estimation by modified IPW estimator:")
      if (any(class(mdl.C) == 'coxph') || any(class(mdl.C) == 'survreg')) {
        message("  strata are based on score from the model of 'mdl.C'.")
        strt = strata.1dim(mdl.C$linear.predictor, m)

      } else if (!is.null(zc.C) && class(zc.C)=='numeric') {
        message("  strata are based on a variable 'zc.C'.")
        strt = strata.1dim(zc.C, m)

      } else if (!is.null(zc.C)) {
        mdl.C = try( {coxph(Surv(time, event==0) ~ ., as.data.frame(zc.C))}, silent=TRUE)
        if (class(mdl.C) == 'try-error')
          message("  fail to obtain Cox model for C. standard IPW estimator is used.")
        if (all(predict(mdl.C) == 0))
          message("  fail to obtain Cox model for C. standard IPW estimator is used.")

        message("  strata are based on score from Cox model by 'zc.C'.")
        strt = strata.1dim(mdl.C$linear.predictor, m)
        ret = c(ret, list(mdl.C=mdl.C))

      } else {
        stop("\n\nnow coxph and survreg object are accepted for censoring time model 'mdl.C'. for other objects, score vecotr for each subject 'zc.C' is to be specified.\n\n")
      }
    } else if (point.method=='ipw') {
      message("\nNRI estimation by standard IPW estimator:")
      strt = rep(1, length(time))
    }
    ret = c(ret, list(strt=strt))
    est = nricens.ipw.main(time, event, upp, dwn, t0, strt)

  } else if (point.method=='km' || point.method=='modkm') {

    if (point.method=='modkm') {
      message("\nNRI estimation by modified KM estimator:")
      message("  strata for all subjects:")
      wk.all = strata.cox.2dim(time, event, zc.T, zc.C, m, pca2nd, !upp & !dwn)
      message("  strata for up subjects:")
      wk.upp = strata.cox.2dim(time, event, zc.T, zc.C, m, pca2nd, upp)
      message("  strata for down subjects:")
      wk.dwn = strata.cox.2dim(time, event, zc.T, zc.C, m, pca2nd, dwn)
      strt     = wk.all$strt
      strt.upp = wk.upp$strt
      strt.dwn = wk.dwn$strt

      ret = c(ret, list(strt.all=strt, strt.up=strt.upp, strt.down=strt.dwn,
        mdl.T.all=wk.all$mdl.T,  mdl.C.all=wk.all$mdl.C,
        mdl.T.up=wk.upp$mdl.T,   mdl.C.up=wk.upp$mdl.C,
        mdl.T.down=wk.dwn$mdl.T, mdl.C.down=wk.dwn$mdl.C,
        pca.all=wk.all$pca, pca.up=wk.upp$pca, pca.down=wk.dwn$pca))
    } else {
      message("\nNRI estimation by standard KM estimator:")
      strt.upp = strt.dwn = strt = rep(1, length(time))
    }
    est = nricens.km.main(time, event, upp, dwn, t0, strt, strt.upp, strt.dwn)
  }
  message("\nPoint estimates:")
  result = data.frame(est)
  names(result) = 'Estimate'
  row.names(result) = c('NRI','NRI+','NRI-','Pr(Up|Case)','Pr(Down|Case)','Pr(Down|Ctrl)','Pr(Up|Ctrl)')
  print(result)

  ##
  ## interval estimation
  if (niter > 0) {
    message("\nNow in bootstrap..")
    N    = length(time)
    samp = matrix(NA, niter, 7)
    colnames(samp) = c('NRI','NRI+','NRI-','Pr(Up|Case)','Pr(Down|Case)','Pr(Down|Ctrl)','Pr(Up|Ctrl)')

    if (point.method == 'ipw' || point.method=='modipw') {
      for (b in 1:niter) {
        f    = as.integer(runif(N, 0, N)) + 1
        objs = list(mdl.std, mdl.new, z.std[f,], z.new[f,], p.std[f], p.new[f])
        wk   = get.uppdwn(time[f], event[f], objs, flag.mdl, flag.prd, flag.rsk, t0, updown, cut, get.risk, msg=FALSE)
        upp  = wk[[1]]
        dwn  = wk[[2]]
        samp[b,] = nricens.ipw.main(time[f], event[f], upp, dwn, t0, strt[f])
      }
    } else if (point.method=='km' || point.method=='modkm') {
      for (b in 1:niter) {
        f    = as.integer(runif(N, 0, N)) + 1
        objs = list(mdl.std, mdl.new, z.std[f,], z.new[f,], p.std[f], p.new[f])
        wk   = get.uppdwn(time[f], event[f], objs, flag.mdl, flag.prd, flag.rsk, t0, updown, cut, get.risk, msg=FALSE)
        upp  = wk[[1]]
        dwn  = wk[[2]]
        samp[b,] = nricens.km.main(time[f], event[f], upp, dwn, t0, strt[f], strt.upp[f], strt.dwn[f])
      }
    }
    ret = c(ret, list(bootstrapsample=samp))

    ci = as.numeric(apply(samp, 2, quantile, c(alpha/2, 1-alpha/2), na.rm=TRUE, type=2))
    message("\nPoint & Interval estimates:")
    result = as.data.frame(cbind(est, matrix(ci, ncol=2, byrow=TRUE)))
    names(result) = c('Estimate', 'Lower', 'Upper')
    row.names(result) = c('NRI','NRI+','NRI-','Pr(Up|Case)','Pr(Down|Case)','Pr(Down|Ctrl)','Pr(Up|Ctrl)')
    print(result)
  }

  invisible(c(list(nri=result), ret))
}
