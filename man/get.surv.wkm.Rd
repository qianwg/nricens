\name{get.surv.wkm}
\alias{get.surv.wkm}
\title{
  Calculate a survival probability at a given time by the weighted
  Kaplan-Meier estimator.
}
\description{
  This function estimates \eqn{Pr\{T > t_0\}} by the weighted
  Kaplan-Meier estimator.
  It is required to specify \eqn{t_0} and strata for the calculation.
  Strata can be obtained by the function of \code{strata.cox.2dim} in
  this package. See examples below.
}
\usage{
get.surv.wkm(time, event, t0, strt)
}
\arguments{
  \item{time}{
    Vector of follow up times.
  }
  \item{event}{
    Vector of event indicators, 1 for event of interest, 0 for censoring.
  }
  \item{t0}{
    Scalar value indicating a time to determine evnet/non-event (\eqn{t_0}).
  }
  \item{strt}{
    Vector of integer values corresponding to the strata number for each
    subject.
  }
}
\examples{
## here consider pbc dataset in survival package as an example
## calculate survival probability at several times by the weighted Kaplan-Meier estimator
library(survival)
dat   = pbc[1:312,]
time  = dat$time
event = ifelse(dat$status==2, 1, 0)
N     = nrow(dat)
niter = 30
alpha = 0.05

## covariates for the event and the censoring time models: age, bilirubin, albumin, and protime
zc = as.matrix(subset(dat, select = c(age, bili, albumin, protime)))

## obtaining strata by the function provided by this package
strt = strata.cox.2dim(time, event, zc, zc, m = 100, pca2nd = 1)$strt

## weighted Kaplan-Meier estimation
results = NULL
ts = seq(min(time), max(time), by=1000)
for (t0 in ts) {
  sv = get.surv.wkm(time, event, t0, strt)

  ## bootstrap
  ## number of iterations should be increased in practice
  samp = rep(NA, niter)
  for (b in 1:niter) {
    f = as.integer(runif(N, 0, N)) + 1
    samp[b] = get.surv.wkm(time[f], event[f], t0, strt[f])
  }
  ci = as.numeric(quantile(samp, c(alpha/2, 1-alpha/2), type=2))

  results = rbind(results, data.frame(t=t0, surv=sv, lower=ci[1], upper=ci[2]))
}
results

## comparison with the standard Kaplan-Merier
summary(survfit(Surv(time, event) ~ 1, dat), t=ts)
}
