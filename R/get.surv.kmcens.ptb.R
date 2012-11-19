get.surv.kmcens.ptb <-
function (kmce, V, t0) {
  N = nrow(kmce$psi)
  return(stepfun(kmce$dtm, c(1, kmce$kap - (t(kmce$psi) %*% (V-1) / N)))(t0))
}
