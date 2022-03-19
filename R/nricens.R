nricens <-
function (time=NULL, event=NULL, mdl.std=NULL, mdl.new=NULL, z.std=NULL, z.new=NULL, p.std=NULL, p.new=NULL,
  t0=NULL, updown='category', cut=NULL, point.method='km', niter=1000, alpha=0.05, msg=TRUE) {

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
  if (point.method=='ipw') {
    message("\nNRI estimation by IPW estimator:")
    est = nricens.ipw.main(time, event, upp, dwn, t0)

  } else if (point.method=='km') {
    message("\nNRI estimation by KM estimator:")
    est = nricens.km.main(time, event, upp, dwn, t0)
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

    if (point.method == 'ipw') {
      for (b in 1:niter) {
        f    = as.integer(runif(N, 0, N)) + 1
        objs = list(mdl.std, mdl.new, z.std[f,], z.new[f,], p.std[f], p.new[f])
        wk   = get.uppdwn(time[f], event[f], objs, flag.mdl, flag.prd, flag.rsk, t0, updown, cut, get.risk, msg=FALSE)
        upp  = wk[[1]]
        dwn  = wk[[2]]
        samp[b,] = nricens.ipw.main(time[f], event[f], upp, dwn, t0)
      }
    } else if (point.method=='km') {
      for (b in 1:niter) {
        f    = as.integer(runif(N, 0, N)) + 1
        objs = list(mdl.std, mdl.new, z.std[f,], z.new[f,], p.std[f], p.new[f])
        wk   = get.uppdwn(time[f], event[f], objs, flag.mdl, flag.prd, flag.rsk, t0, updown, cut, get.risk, msg=FALSE)
        upp  = wk[[1]]
        dwn  = wk[[2]]
        samp[b,] = nricens.km.main(time[f], event[f], upp, dwn, t0)
      }
    }
    ret = c(ret, list(bootstrapsample=samp))

    ci = as.numeric(apply(samp, 2, quantile, c(alpha/2, 1-alpha/2), na.rm=TRUE, type=2))
    se = as.numeric(apply(samp, 2, sd))
    pvalue = (1-pnorm(abs(est/se))) *2
    message("\nPoint & Interval estimates:")
    result = as.data.frame(cbind(est, se, matrix(ci, ncol=2, byrow=TRUE), pvalue))
    names(result) = c('Estimate', 'Std.Error', 'Lower', 'Upper', 'Pvalue')
    row.names(result) = c('NRI','NRI+','NRI-','Pr(Up|Case)','Pr(Down|Case)','Pr(Down|Ctrl)','Pr(Up|Ctrl)')
    print(result)
  }

  invisible(c(list(nri=result), ret))
}
