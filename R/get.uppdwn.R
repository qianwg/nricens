get.uppdwn <-
function (time, event, z.std, z.new, t0, cut, region.method, msg=F) {

  mdl.std = coxph(Surv(time, event==1) ~ ., as.data.frame(z.std))
  mdl.new = coxph(Surv(time, event==1) ~ ., as.data.frame(z.new))

  p.std = get.risk.cox(mdl.std, t0)
  p.new = get.risk.cox(mdl.new, t0)

  if (region.method == 'cutoff') {
    p.stdc = categorize(p.std, threshold=cut)
    p.newc = categorize(p.new, threshold=cut)
    upp = p.newc - p.stdc > 0
    dwn = p.newc - p.stdc < 0
  } else if (region.method == 'diff') {
    upp = ifelse(p.new - p.std > cut, T, F)
    dwn = ifelse(p.std - p.new > cut, T, F)
  }

  ## message
  if (msg) {
    case = time <= t0 & event == 1
    ctrl = time >  t0

    if (region.method == 'cutoff') {
      print(table(p.stdc, p.newc))
      print(table(p.stdc[case], p.newc[case]))
      print(table(p.stdc[ctrl], p.newc[ctrl]))
    } else if (region.method == 'diff') {
      print(table(upp));       print(table(dwn))
      print(table(upp[case])); print(table(dwn[case]))
      print(table(upp[ctrl])); print(table(dwn[ctrl]))
    }

    plot(p.std, p.new, xlab='Standard', ylab='New', main='')
    lines(0:1000/1000, 0:1000/1000, col=2)
    if (region.method == 'diff') {
      lines(0:1000/1000, 0:1000/1000-cut, lty=2, col=4)
      lines(0:1000/1000, 0:1000/1000+cut, lty=2, col=4)
    }

    yl = range(c(p.std, p.new))

    prb = rbind(p.std[case], p.new[case])
    matplot(prb, type='l', xlab='', ylab='', main='Case', axes=F, xlim=c(0.85,2.15), ylim=yl)
    axis(1, at=c(1,2), labels=c('Standard','New')); axis(2); box(lwd=2)

    prb = rbind(p.std[ctrl], p.new[ctrl])
    matplot(prb, type='l', xlab='', ylab='', main='Control', axes=F, xlim=c(0.85,2.15), ylim=yl)
    axis(1, at=c(1,2), labels=c('Standard','New')); axis(2); box(lwd=2)
  }

  return(list(upp, dwn))
}
