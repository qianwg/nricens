\name{get.surv.wkm}
\alias{get.surv.wkm}
\title{
  Calculating the survival probability by the weighted Kaplan-Meier
  estimator at the given time.
}
\description{
  Calculating the survival probability by the standard Kaplan-Meier
  estimator given the time. Subgroup definition should be given for calculation.
}
\usage{
get.surv.wkm(time, event, t0, strt)
}
\arguments{
  \item{time}{
    Vector of time to event
  }
  \item{event}{
    Vector of event indicator, 1 for event of interest, 0 for censoring
  }
  \item{t0}{
    Time of defining evnet/non-event
  }
  \item{strt}{
    Vector values of the subgroup definition for weighted Kaplan-Meier estimator
  }
}
\examples{
##
## an example to calculate survival probability by the weighted Kaplan-Meier estimator

library(survival)

dat   = na.omit(pbc[1:312,])
time  = dat$time
event = as.integer(dat$status == 2)
N     = nrow(dat)
zcT   = subset(dat, select=c(age, sex, albumin, bili, protime))
zcC   = subset(dat, select=c(age, sex, albumin, bili, protime))
niter = 50
alpha = 0.05

strt = strata.cox.2dim(time, event, zcT, zcC, ms.strt = 100, pca2nd = 1)

## weighted Kaplan-Meier estimation
results = NULL
for (t0 in seq(min(time), max(time), by=1000)) {
  sv = get.surv.wkm(time, event, t0, strt)

  ## bootstrap
  samp = rep(NA, niter)
  for (b in 1:niter) {
    f = as.integer(runif(N, 0, N)) + 1
    samp[b] = get.surv.wkm(time[f], event[f], t0, strt[f])
  }
  ci = as.numeric(quantile(samp, c(alpha/2, 1-alpha/2), type=2))

  results = rbind(results, data.frame(t=t0, surv=sv, lower=ci[1], upper=ci[2]))
}
results
}
