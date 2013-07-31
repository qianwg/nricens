get.uppdwn.bin <-
function (event, objs, flag.mdl, flag.prd, flag.rsk, updown, cut, link, msg=FALSE) {

  mdl.std = objs[[1]]
  mdl.new = objs[[2]]
  z.std   = objs[[3]]
  z.new   = objs[[4]]
  p.std   = objs[[5]]
  p.new   = objs[[6]]

  if (flag.mdl || flag.prd) {
    p.std = update(mdl.std)$fitted.values
    p.new = update(mdl.new)$fitted.values
  }

  rtab = rtab.case = rtab.ctrl = NULL
  if (updown == 'category') {
    p.stdc = categorize(p.std, threshold = cut)
    p.newc = categorize(p.new, threshold = cut)
    upp = p.newc - p.stdc > 0
    dwn = p.newc - p.stdc < 0
  } else if (updown == 'diff') {
    upp = ifelse(p.new - p.std > cut, TRUE, FALSE)
    dwn = ifelse(p.std - p.new > cut, TRUE, FALSE)
  }

  ## message
  if (msg) {
    case = event == 1
    ctrl = event == 0
    message("\nUP and DOWN calculation:")
    message(paste("  #of total, case, and control subjects at t0: ", length(event), sum(case), sum(ctrl)))

    if (updown == 'category') {
      p.stdc    = factor(p.stdc, levels=1:(length(cut)+1), labels=c(paste('<', cut), paste('>=', cut[length(cut)])))
      p.newc    = factor(p.newc, levels=1:(length(cut)+1), labels=c(paste('<', cut), paste('>=', cut[length(cut)])))
      rtab      = table(Standard = p.stdc,       New = p.newc)
      rtab.case = table(Standard = p.stdc[case], New = p.newc[case])
      rtab.ctrl = table(Standard = p.stdc[ctrl], New = p.newc[ctrl])

      message("\n  Reclassification Table for all subjects:")
      print(rtab)
      message("\n  Reclassification Table for case:")
      print(rtab.case)
      message("\n  Reclassification Table for control:")
      print(rtab.ctrl)

    } else if (updown == 'diff') {
      message(paste("  #of subjects with 'p.new - p.std > cut' for all, case, control:", sum(upp), sum(upp[case]), sum(upp[ctrl])))
      message(paste("  #of subjects with 'p.std - p.new < cut' for all, case, control:", sum(dwn), sum(dwn[case]), sum(dwn[ctrl])))
    }

    plot(p.std[case], p.new[case], xlab='Standard model', ylab='New model', main='', xlim=c(0,1), ylim=c(0,1), col=2)
    par(new=T)
    plot(p.std[ctrl], p.new[ctrl], xlab='', ylab='', axes=F, xlim=c(0,1), ylim=c(0,1))
    legend('topleft', c('Case','Control'), col=c(2,1), bty='n', pch=1)
    
    if (updown == 'diff') {
      abline(0,1)
      abline(-cut, 1, lty=2)
      abline(cut, 1, lty=2)
    }
    if (updown == 'category') {
      abline(h=cut, lty=2)
      abline(v=cut, lty=2)
    }
  }

  return(list(upp, dwn, p.std, p.new, rtab, rtab.case, rtab.ctrl))
}
