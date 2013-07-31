\name{nricens.km.main}
\alias{nricens.km.main}
\title{
  Estimate NRI by the standard Kaplan-Meier method.
}
\description{
  Internaly used function by \code{nricens} to provide the NRI estimator
  by the Kaplan-Meier(KM) method.
  This function also provides the NRI estimator by the weighted KM
  method.
  In the comparison for time to event models, it is possible to use
  directly when \code{UP} and \code{DOWN} subjects are known.
}
\usage{
nricens.km.main(time, event, upp, dwn, t0, strt, strt.upp, strt.dwn)
}
\arguments{
  \item{time}{
    Vector of follow up times.
  }
  \item{event}{
    Vector of event indicators, 1 for event of interest, 0 for censoring.
  }
  \item{upp, dwn}{
    Vectors of logical values to determine subjects who belong to
    \code{UP} and \code{DOWN}, respectively.
  }
  \item{t0}{
    Scalar value indicating a time to determine evnet/non-event.
  }
  \item{strt, strt.upp, strt.dwn}{
    Vectors of strata number for each subject.
    These items are used for the estimation by the modified KM NRI
    estimator.
    \code{strt} determines strata for estimating survival function
    for subjects who do not belong to both UP and DOWN,
    \code{strt.upp} determines strata for estimating survival function
    for subjects who belong to UP, and
    \code{strt.dwn} determines strata for estimating survival function
    for subjects who belong to DOWN.
  }
}
