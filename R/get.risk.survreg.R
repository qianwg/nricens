get.risk.survreg <-
function(mdl, t0) {
  dens = survreg.distributions[[survreg.distributions[[mdl$dist]]$dist]]$density
  risk = dens( (-mdl$linear.predictors + log(t0) ) / mdl$scale)[,1]
  
  return(risk)
}
