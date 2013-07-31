\name{nricens.ipw.main}
\alias{nricens.ipw.main}
\title{
  Estimate NRI by the inverse probability weighting (IPW) method.
}
\description{
  Internaly used function by \code{nricens} to provide the NRI estimator
  by the IPW method.
  This function also provides the modified IPW NRI estimator.
  In the comparison for time to event models, it is possible to use
  directly when \code{UP} and \code{DOWN} subjects are known.
}
\usage{
nricens.ipw.main(time, event, upp, dwn, t0, strt)
}
\arguments{
  \item{time}{
    Vector of follow up times.
  }
  \item{event}{
    Vector of event indicators, 1 for event of interest, 0 for censoring.
  }
  \item{upp, dwn}{
    Vector of logical values to determine subjects who belong to
    \code{UP} and \code{DOWN}, respectively.
  }
  \item{t0}{
    Scalar value indicating a time to determine evnet/non-event.
  }
  \item{strt}{
    Vector of strata number for each subject.
    This item is used for the estimation by the modified IPW NRI estimator.
  }
}
