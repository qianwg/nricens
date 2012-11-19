get.risk.cox <-
function(mdl, t0) {
  bash    = basehaz(mdl)
  lambda0 = approx(bash$time, bash$hazard, t0)$y
  risk    = 1 - exp( - lambda0 * exp( mdl$linear.predictors ) )
  return(risk)
}
