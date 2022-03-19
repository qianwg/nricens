nribin <-
function (event=NULL, mdl.std=NULL, mdl.new=NULL, z.std=NULL, z.new=NULL, p.std=NULL, p.new=NULL,
  updown='category', cut=NULL, link='logit', niter=1000, alpha=0.05, msg=TRUE) {

  ##
  ## type of calculation
  flag.mdl = !is.null(mdl.std) && !is.null(mdl.new)
  flag.prd = !is.null(z.std)   && !is.null(z.new)
  flag.rsk = !is.null(p.std)   && !is.null(p.new)

  ##
  ## check standard & new model
  if (flag.mdl) {
    if (is.null(event)) event = as.numeric(mdl.std$y)
    if (is.null(mdl.std$x) || is.null(mdl.new$x))
      stop("\n\nmodel object does not contain predictors. pls set x=TRUE for model calculation.\n\n")

    z.std   = mdl.std$x[,-1]
    z.new   = mdl.new$x[,-1]
    link    = mdl.std$family[[2]]
    mdl.std = glm(event ~ ., family=binomial(link), data=as.data.frame(cbind(event, z.std)))
    mdl.new = glm(event ~ ., family=binomial(link), data=as.data.frame(cbind(event, z.new)))

  } else if (flag.prd) {
    mdl.std = glm(event ~ ., family=binomial(link), data=as.data.frame(cbind(event, z.std)))
    mdl.new = glm(event ~ ., family=binomial(link), data=as.data.frame(cbind(event, z.new)))
    message("\nSTANDARD prediction model:")
    print(summary(mdl.std)$coef)
    message("\nNEW prediction model:")
    print(summary(mdl.new)$coef)

  } else if (!flag.mdl && !flag.prd && !flag.rsk) {
    stop("\n\neither one of 'event, z.std, z.new', 'event, p.std, p.new', and 'mdl.std, mdl.new' should be specified.\n\n")
  }

  if (is.null(cut))
    stop("\n\n'cut' is empty")

  objs = list(mdl.std, mdl.new, z.std, z.new, p.std, p.new)

  ##
  ## DH & DL
  wk  = get.uppdwn.bin(event, objs, flag.mdl, flag.prd, flag.rsk, updown, cut, link, msg=msg)
  upp = wk[[1]]
  dwn = wk[[2]]
  ret = list(mdl.std=mdl.std, mdl.new=mdl.new, p.std=wk[[3]], p.new=wk[[4]], up=upp, down=dwn, rtab=wk[[5]], rtab.case=wk[[6]], rtab.ctrl=wk[[7]])

  ##
  ## point estimation
  message("\nNRI estimation:")
  est = nribin.count.main(event, upp, dwn)
  message("Point estimates:")
  result = data.frame(est)
  names(result) = 'Estimate'
  row.names(result) = c('NRI','NRI+','NRI-','Pr(Up|Case)','Pr(Down|Case)','Pr(Down|Ctrl)','Pr(Up|Ctrl)')
  print(result)

  ##
  ## interval estimation
  if (niter > 0) {
    message("\nNow in bootstrap..")
    ci   = rep(NA, 14)
    N    = length(event)
    samp = matrix(NA, niter, 7)
    colnames(samp) = c('NRI','NRI+','NRI-','Pr(Up|Case)','Pr(Down|Case)','Pr(Down|Ctrl)','Pr(Up|Ctrl)')

    for (b in 1:niter) {
      f    = as.integer(runif(N, 0, N)) + 1
      objs = list(mdl.std, mdl.new, z.std[f,], z.new[f,], p.std[f], p.new[f])
      wk   = get.uppdwn.bin(event[f], objs, flag.mdl, flag.prd, flag.rsk, updown, cut, link, msg=FALSE)
      upp  = wk[[1]]
      dwn  = wk[[2]]
      samp[b,] = nribin.count.main(event[f], upp, dwn)
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
