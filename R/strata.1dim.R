strata.1dim <-
function (score, m, subs=NULL) {

  n    = length(score)
  strt = rep(1, n)

  if (is.null(subs)) subs = rep(TRUE, n)

  score.w = score[subs]
  nr      = length(score.w)
  nstr    = as.integer(nr / m)

  if (is.null(score) || nstr <= 1) return(strt)

  strt = categorize(score, threshold = quantile(score.w, (1:nstr)/nstr))

  return(strt)
}
