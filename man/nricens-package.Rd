\name{nricens-package}
\alias{nricens-package}
\docType{package}
\title{
  NRI for risk prediction models with time to event data
}
\description{
  Calculating NRI for comparing risk prediction models for time to event
  outcome. The Cox regression model is used as a prediction model.
  The modified estimators to deal with an independent censoring assumption
  are provided. Confidence interval for each NRI term are provided by the
  percentile bootstrap method.
}
\details{
\tabular{ll}{
Package: \tab nricens\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-11-19\cr
License: \tab GPL-2\cr
}
}
\author{
Eisuke Inoue

Maintainer: Eisuke Inoue <inouee@pharm.kitasato-u.ac.jp>
}
\examples{
library(survival)

dat = na.omit(pbc[1:312,])
dat$sex = ifelse(dat$sex=='f', 1, 0)

## time to event
time = dat$time

## event indicator
event = as.integer(dat$status == 2)

## predictors in the standard model
z.std = subset(dat, select=c(age, sex, bili, chol))

## predictors in the new model
z.new = subset(dat, select=c(age, sex, bili, chol, protime))

## time to define 'event' / 'non-event'
t0 = 2000

## for standard Kaplan-Meier based NRI estimator
#nricens(time, event, z.std, z.new, t0, cut=0.05, niter=1000, region.method='diff')
nricens(time, event, z.std, z.new, t0, cut=0.05, niter=50, region.method='diff')

## for modified Kaplan-Meier based NRI estimator
#nricens(time, event, z.std, z.new, t0, cut=0.05, niter=1000, zcT=z.new, ms.strt=100, region.method='diff')

## for standard IPW based NRI estimator
#nricens(time, event, z.std, z.new, t0, cut=0.05, niter=1000, region.method='diff', point.method='ipw')

## for modified Kaplan-Meier based NRI estimator
#nricens(time, event, z.std, z.new, t0, cut=0.05, niter=1000, zcC=z.new, ms.strt=100, region.method='diff', point.method='ipw')

}
