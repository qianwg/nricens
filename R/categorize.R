categorize <-
function(dat, threshold=NULL, nlevel=NULL) {

  if (is.null(threshold))
    threshold = quantile(dat, seq(0,1,1/nlevel), na.rm=T)[-1]
  else
    nlevel = length(threshold)

  threshold[nlevel] = max(dat, na.rm=T) + 1

  new = rep(NA, length(dat))
  for (i in 1:nlevel)
    new[is.na(new) & dat < threshold[i]] = i

  return(new)
}
