\name{nricens-package}
\alias{nricens-package}
\alias{nricens-package}
\docType{package}
\title{
  R Functions to calculate NRI for comparing time to event and binary
  response models.
}
\description{
  This package provides the functions to estimate the net
  reclassification improvement (NRI) for competing risk prediction
  models with time to event and binary response data.
  The NRI for binary response models can be calculated by \code{nribin},
  and that for time to event models can be calculated by
  \code{nricens}.
  The risk category based NRI and the risk difference based NRI are
  provided by these functions.
  Users can use several estimators for comparing time to event models.
  Confidence intervals are calculated by the percentile bootstrap
  method.
  In this version, several types of input data are allowed to improve
  user-friendliness and convenience.
}
\details{
\tabular{ll}{
Package: \tab nricens\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2013-07-31\cr
License: \tab GPL-2\cr
}
}
\author{
Eisuke Inoue <inouee@pharm.kitasato-u.ac.jp>
}
\references{
Pencina MJ, D'Agostino RB, Steyerberg EW. Extensions of net
reclassification improvement calculations to measure usefulness of new
biomarkers. Statistics in Medicine 2011.

Uno H, Tian L, Cai T, Kohane IS, Wei LJ. A unified inference procedure
for a class of measures to assess improvement in risk prediction systems
with survival data, Statistics in Medicine 2012.

Hsu CH, Taylor JMG. A robust weighted Kaplan-Meier approach for data
with dependent censoring using linear combinations of prognostic
covariates, Statistics in Medicine 2010.
}
\examples{
## consider pbc dataset in survival package as an example
library(survival)
dat = pbc[1:312,]
dat$sex = ifelse(dat$sex=='f', 1, 0)

## predciting the event of 'death'
time  = dat$time
event = ifelse(dat$status==2, 1, 0)

## standard prediction model: age, bilirubin, and albumin
z.std = as.matrix(subset(dat, select = c(age, bili, albumin)))

## new prediction model: age, bilirubin, albumin, and protime
z.new = as.matrix(subset(dat, select = c(age, bili, albumin, protime)))

## fitting cox models
mstd = coxph(Surv(time, event) ~ ., data.frame(time, event, z.std), x=TRUE)
mnew = coxph(Surv(time, event) ~ ., data.frame(time, event, z.new), x=TRUE)

## Calculation of the risk category NRI at 2000 days
nricens(mdl.std = mstd, mdl.new = mnew, t0 = 2000, cut = c(0.2, 0.4),
        niter = 10)


## Next, consider binary prediction models
library(survival)
dat = pbc[1:312,]
dat$sex = ifelse(dat$sex=='f', 1, 0)

## subjects censored before 2000 days are excluded
dat = dat[ dat$time > 2000 | (dat$time < 2000 & dat$status == 2), ]

## predciting the event of 'death' before 2000 days
event = ifelse(dat$time < 2000 & dat$status == 2, 1, 0)

## standard prediction model: age, bilirubin, and albumin
z.std = as.matrix(subset(dat, select = c(age, bili, albumin)))

## new prediction model: age, bilirubin, albumin, and protime
z.new = as.matrix(subset(dat, select = c(age, bili, albumin, protime)))

## glm fit (logistic model)
mstd = glm(event ~ ., binomial(logit), data.frame(event, z.std), x=TRUE)
mnew = glm(event ~ ., binomial(logit), data.frame(event, z.new), x=TRUE)

## Calculation of risk difference NRI
nribin(mdl.std = mstd, mdl.new = mnew, cut = 0.02, niter = 0,
       updown = 'diff')
}
