\name{nricens}
\alias{nricens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function of the \code{nricens} package
}
\description{
Calculating NRI for comparing risk prediction models for time to event
outcome. The Cox regression model is used as a prediction model.
The modified estimators to deal with an independent censoring assumption
are provided. Confidence interval for each NRI term are provided by the
percentile bootstrap method.
}
\usage{
nricens(time, event, z.std, z.new, t0, cut, niter = 1000, zcT = NULL,
zcC = NULL, ms.strt = 30, pca2nd = 1, region.method = "cutoff",
point.method = "km",  alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
    Vector of time to event
}
  \item{event}{
    Vector of event indicator, 1 for event of interest, 0 for censoring
}
  \item{z.std}{
    Matrix of predictors for a standard model. Factor or character
    variables are not allowed. Missing values are not allowed.
}
  \item{z.new}{
    Matrix of predictors for new or improved model.
}
  \item{t0}{
    Time of defining evnet/non-event.
}
  \item{cut}{
    Scalar or vector value.
    For a risk-cutoff NRI calculation, \code{cut} is a scalar or vector
    value to define risk cutoff value. For a risk-difference NRI
    calculation, \code{cut} is a scalar to define the threshold of the
    difference in risk. When \code{cut=0} is specified, continuous NRI
    will be provided. See also \code{region.method}.
}
  \item{niter}{
    Scalar value to define the number of iterations for bootstratp
    calculation.
    When \code{niter=0}, calculation for a confidence interval is skipped.
}
  \item{zcT}{
    Matrix of covariates for a working time-to-event model which is used
    for modified Kaplan-Meier-based NRI estimator.
    The default value is \code{NULL}, which will provide the estimate by
    the standard Kaplan-Meier-based NRI estimator.
    When it is not NULL, modified NRI estimator will be applied.
}
  \item{zcC}{
    Matrix of covariates for a working censoring time model which is
    used for modified KM-based and modified IPW-based NRI estimators.
    The default value is \code{NULL}.
    Wnen \code{point.method='km'}, \code{zcT} is not \code{NULL}, and
    \code{zcC = NULL}, the same matrix of covariates of \code{zcT} is
    used for \code{zcC}.
}
  \item{ms.strt}{
    Scalar to define the number of minimum samples in a strata.
}
  \item{pca2nd}{
    Scalar for the number of strata of the second principle component.
}
  \item{region.method}{
    Character to define the region for detemining \code{UP} and
    \code{DOWN} in NRI calculation. For \code{"region.method = 'cutoff'"},
    the risk cutoff based NRI is calculated. For \code{"region.method = 'diff'"},
    the risk difference based NRI is calculated.
}
  \item{point.method}{
    Character to define the type of estimator.
    For \code{"point.method = 'km'"}, the standard Kaplan-Meier based
    estimator will be applied.
    For \code{"point.method = 'ipw'"}, the inverse probability for
    censoring weighted estimator will be applied.
}
  \item{alpha}{
    1-alpha confidence interval is calcualted.
}
}
