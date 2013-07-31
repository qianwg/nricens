categorize <-
function(dat, threshold) {

  threshold[length(threshold) + 1] = max(dat, na.rm=TRUE) + 1

  new = rep(NA, length(dat))
  for (i in 1:length(threshold))
    new[is.na(new) & dat < threshold[i]] = i

  return(new)
}
