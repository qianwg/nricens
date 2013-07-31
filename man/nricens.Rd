\name{nricens}
\alias{nricens}
\title{
NRI for time to event models
}
\description{
  This function estimates the NRI for competing risk prediction models
  with time to event variable. \code{coxph} object, \code{survreg}
  object, predictors, and predicted risks can be used as input data for
  the calculation.
  The risk category based NRI and the risk difference based NRI can be
  calculated.
  Users can use several types of estimators to obtain point estimates of
  the NRI and its components.
  The percentile bootstrap method is used for an interval estimation.
}
\usage{
nricens(time = NULL, event = NULL, mdl.std = NULL, mdl.new = NULL,
        z.std = NULL, z.new = NULL, p.std = NULL, p.new = NULL, t0 = NULL,
        updown = "category", cut = NULL, point.method = "km",
        zc.T = NULL, zc.C = NULL, mdl.C = NULL, m = 100, pca2nd = 1,
        niter = 1000, alpha = 0.05, msg = TRUE)
}
\arguments{
  \item{time}{
    Vector of observed follow up times, \eqn{X = \min(T, C)}.
    \eqn{T} is event time, and \eqn{C} is censoring time.
  }
  \item{event}{
    Vector of event indicators, 1 for event of interest, 0 for censoring.
  }
  \item{mdl.std, mdl.new}{
    \code{coxph} or \code{survreg} objects corresponding to a standard
    and a new risk prediction model, respectively.
    Since predictors are extracted from these objects, \code{x=TRUE} is
    required when fitting a Cox or a parametric survival model.
  }
  \item{z.std, z.new}{
    Matrix of predictors for a standard and a new risk prediction
    model, respectively.
    Neither factor nor character nor missing values are allowed.
  }
  \item{p.std, p.new}{
    Vector of predicted risks from a standard and a new prediction
    model, respectively.
  }
  \item{t0}{
    Scalar value indicating a time to determine evnet/non-event.
  }
  \item{updown}{
    Character to specify the method to determine \code{UP} and
    \code{DOWN}.
    When \code{"category"} is specified (by default), the risk category
    based NRI is calculated.
    When \code{"diff"} is specified, the risk difference based NRI is
    calculated.
  }
  \item{cut}{
    Scalar or vector values to specify the cutoff value(s) of predicted
    risks for determining \code{UP} and \code{DOWN}.
    For the risk category based NRI (\code{updown = "category"}), this
    option corresponds to cutoff value(s) of risk categories.
    For the risk difference based NRI (\code{updown = "diff"}), this
    option corresponds to a cutoff value of a risk difference
    (only scalar is allowed).
  }
  \item{point.method}{
    Character to determine an estimator for a point estimation.
    When \code{"km"} is specified, the Kaplan-Meier (KM) based NRI
    estimator is used (Pencina et al., 2011). 
    When \code{"ipw"} is specified, the inverse probability weighting
    (IPW) NRI estimator is used (Uno et al., 2012).
    When \code{"modkm"} is specified, the modified KM based NRI
    estimator is used.
    When \code{"modipw"} is specified, the modified IPW NRI estimator is
    used.
  }
  \item{zc.T}{
    Matrix of covariates for a working Cox model for \eqn{T},
    which is used to construct strata for the modified KM based NRI
    estimator.
  }
  \item{zc.C}{
    Matrix of covariates for a working Cox model for \eqn{C},
    which is used to construct strata for the modified KM and the
    modified IPW NRI estimators.
  }
  \item{mdl.C}{
    \code{coxph} or \code{survreg} object to construct strata for the
    modified IPW NRI estimator.
  }
  \item{m}{
    Scalar value to determine the number of minimum samples for each
    stratum, which is used both in the modified KM and the modified IPW
    NRI estimators.
  }
  \item{pca2nd}{
    Scalar value to determine the number of strata in the second
    principal component, which is used in the modified KM based NRI
    estimator.
  }
  \item{niter}{
    Scalar value to determine the number of bootstratp sampling. When
    \code{0} is specified, an interval estimation is skipped.
  }
  \item{alpha}{
    1-alpha confidence interval is calcualted.
  }
  \item{msg}{
    Logical value to display computation process.
    Setting \code{FALSE} leads silent execution.
  }
}
\details{
  Either one set of the following arguments should be specified for the
  NRI calculation: (\code{mdl.std}, \code{mdl.new});
  (\code{time}, \code{event}, \code{z.std}, \code{z.new});
  and (\code{time}, \code{event}, \code{p.std}, \code{p.new}).

  In the first set of the argument, (\code{mdl.std}, \code{mdl.new}),
  fitted results by \code{coxph} or \code{survreg} are used for the NRI
  calculation.
  \code{time}, \code{event}, \code{z.std}, and \code{z.new} are
  extracted from fitted result objects.
  The variance of model parameters are accounted for an interval
  estimation of the NRI.
  When \code{time} and \code{event} are specified in arguments,
  those specified are used without extracting from \code{coxph} or
  \code{survreg} objects.

  In the second set of the argument, (\code{time}, \code{event},
  \code{z.std}, \code{z.new}), a standard and a new prediction models
  are fitted inside this function with \code{time}, \code{event},
  \code{z.std} and \code{z.new}.
  The variance of model parameters are also accounted for an interval
  estimation of the NRI.

  In the third set of the argument, (\code{time}, \code{event},
  \code{p.std}, \code{p.new}), predicted risks are used.
  Since fit of prediction models are not conducted while in a bootstrap,
  this can be used for a validation study by an external data source or
  by a cross-validation.

  For the risk category based NRI calculation, cutoff values of risk
  category can be specified by \code{cut}, which is a scalar for the
  case of two risk categories and is a vector for the case of more than
  two risk categories.
  \code{UP} and \code{DOWN} are determined by the movement in risk
  categories.

  For the risk difference based NRI calculation, cutoff values of risk
  difference can also be specified by \code{cut}, where \code{UP} and
  \code{DOWN} are defiend as \eqn{p_{new} - p_{standard} > \delta} and
  \eqn{p_{standard} - p_{new} > \delta}, respectively.
  \eqn{p_{standard}} and \eqn{p_{new}} are predicted individual risks
  from a standard and a new prediction model, respectively, and
  \eqn{\delta} corresponds to \code{cut}.
  The continuous NRI, which is the special version of the risk
  difference based NRI, can be calculated by specifying both
  \code{updown = "diff"} and \code{cut = 0}.

  When the modified KM based NRI estimator
  (\code{point.method = "modkm"}) is specified, the weighted KM
  estimators (Hsu and Taylor, 2010) instead of the standard KM are used
  for estimating the survival functions in the NRI prospective form
  (Pencina et al., 2011).
  Strata are constructed based on a categorization of individual scores
  from working Cox models for \eqn{T} and \eqn{C}.
  Principal component analyses are conducted before a categorization
  (Hsu and Taylor, 2010).

  When the modified IPW NRI estimator (\code{point.method = "modipw"})
  is specified, the stratified inverse censoring probability weight is
  used.
  Strata are constructed based on a categorization of individual scores
  from a working censoring time model.
  Weights are calculated by the standard Kaplan-Meier estimator
  separately for each stratum.
  Note that weights are not the inverse of predicted censoring time from
  a working model.
  Fitted result of the \code{coxph} or \code{survreg} for \eqn{C},
  covariates of a working Cox model for \eqn{C}, and scores derived from
  other regression models for \eqn{C} can be specified for a
  stratification.
  When covariates of a working model for \eqn{C} are specified to
  \code{zc.C}, a Cox model is fitted inside this function.
  When one dimensional scores (numeric vector) are specified to
  \code{zc.C}, just a categorization is conducted to create strata.

  Interval estimation is based on the percentile bootstrap method.
}
\value{
  Returns a list of the following items:
  \item{nri}{
    Point and interval estimates of the NRI and its components.
  }
  \item{mdl.std, mdl.new}{
    Fitted \code{coxph} or \code{survreg} objects corresponding to a
    standard and a new prediction model, respectively.
    These items are provided when prediction models are fitted inside
    this function. Otherwise \code{NULL}.
  }
  \item{z.std, z.new}{
    Predictors of a standard and a new prediction model, respectively.
    These items are provided when they are extracted from \code{mdl.std}
    and \code{mdl.new}. Otherwise \code{NULL}.
  }
  \item{p.std, p.new}{
    Predicted risks by a standard and a new prediction model,
    respectively.
  }
  \item{up, down}{
    Logical values to show subjects who belong to \code{UP} and
    \code{DOWN}, respectively.
  }
  \item{rtab, rtab.case, rtab.ctrl}{
    \code{table} objects corresponding to reclassification tables for all,
    case, and control subjects, respectively.
    These items are provided only when the risk category based NRI is
    specified and \code{msg = TRUE}.
  }
  \item{strt.all, strt.up, strt.down}{
    Vectors of strata numbers for each subject to estimate survival
    functions by the weighted Kaplan-Meier method.
    These items are presented only when the modified KM based NRI
    estimator is applied.
    \code{strt.all} shows a strata number to estimate a marginal
    survival function,
    \code{strt.up} shows a strata number to estimate a survival function
    for \code{UP}, and
    \code{strt.down} shows a strata number to estimate a survival
    function for \code{DOWN}.
  }
  \item{mdl.T.all, mdl.C.all, mdl.T.up, mdl.C.up, mdl.T.down, mdl.C.down}{
    Fitted results of the working Cox models for \eqn{T} and \eqn{C}.
    These items are presented only when the modified KM based NRI
    estimator is applied.
    \code{mdl.*.all} is fitted results of the Cox model by subjects who do
    not belong to both \code{UP} and \code{DOWN},
    \code{mdl.*.up} is fitted results of the Cox model by subjects who
    belong to \code{UP}, and
    \code{mdl.*.down} is fitted results of the Cox model by subjects who
    belong to \code{DOWN}.
  }
  \item{pca.all, pca.up, pca.down}{
    Results of the principal component analyses by \code{prcomp}
    using scores from working Cox models for \eqn{T} and \eqn{C}.
    These items are presented only when the modified KM based NRI
    estimator is applied.
    \code{pca.all} is a result of \code{prcomp} by subjects who do not
    belong to both \code{UP} and \code{DOWN},
    \code{pca.up} is a result of \code{prcomp} by subjects who belong to
    \code{UP}, and
    \code{pca.down} is result of \code{prcomp} by subjects who belong to
    \code{DOWN}.
  }
  \item{mdl.C}{
    Fitted result of a working Cox model for \eqn{C} with predictors
    specified to \code{zc.C}.
    This item is presented only when the modified IPW NRI estimator is
    applied.
  }
  \item{strt}{
    Vector of strata number for each subject.
    This item is presented only when the modified IPW NRI estimator is
    applied.
  }
  \item{bootstrapsample}{
    Results of each bootstrap sample.
  }
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
## here consider pbc dataset in survival package as an example
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

## coxph fit
mstd = coxph(Surv(time,event) ~ ., data.frame(time,event,z.std), x=TRUE)
mnew = coxph(Surv(time,event) ~ ., data.frame(time,event,z.new), x=TRUE)

## predicted risk at t0=2000
p.std = get.risk.coxph(mstd, t0=2000)
p.new = get.risk.coxph(mnew, t0=2000)

## Calculation of risk category NRI
##   by the KM estimator using ('mdl.std', 'mdl.std').
nricens(mdl.std = mstd, mdl.new = mnew, t0 = 2000, cut = c(0.2, 0.4),
        niter = 0)

##   by the KM estimator using ('time', 'event', 'z.std', 'z.std').
nricens(time = time, event = event, z.std = z.std, z.new = z.new,
        t0 = 2000, cut = c(0.2, 0.4), niter = 0)

##   by the KM estimator using ('time','event','p.std','p.std').
nricens(time = time, event = event, p.std = p.std, p.new = p.new,
        t0 = 2000, cut = c(0.2, 0.4), niter = 0)

## Calculation of risk difference NRI by the KM estimator
nricens(mdl.std = mstd, mdl.new = mnew, t0 = 2000, updown = 'diff',
        cut = 0.05, niter = 0)

## Calculation of risk difference NRI by the IPW estimator
nricens(mdl.std = mstd, mdl.new = mnew, t0 = 2000, updown = 'diff',
        cut = 0.05, point.method = 'ipw', niter = 0)

## Calculation of risk difference NRI by the modified KM estimator
nricens(mdl.std = mstd, mdl.new = mnew, t0 = 2000, updown = 'diff',
        cut = 0.05, point.method = 'modkm', zc.T = z.new, zc.C = z.new,
        m = 20, niter = 0)

## Calculation of risk difference NRI by the modified IPW estimator
nricens(mdl.std = mstd, mdl.new = mnew, t0 = 2000, updown = 'diff',
        cut = 0.05, point.method = 'modipw', zc.C = z.new,
        m = 100, niter = 0)

## Calculation of risk difference NRI by the modified IPW estimator
## censoring time model is used for constructing strata
mdl.C = coxph(Surv(time, event==0) ~ ., data.frame(time, event, z.new))
nricens(mdl.std = mstd, mdl.new = mnew, t0 = 2000, updown = 'diff',
        cut = 0.05, point.method = 'modipw', mdl.C = mdl.C,
        m = 100, niter = 0)
}
