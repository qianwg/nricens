\name{nricens.ipw.main}
\alias{nricens.ipw.main}
\title{
  Estimate NRI by the inverse probability weighting (IPW) method.
}
\description{
  Internaly used function by \code{nricens} to provide the NRI estimator
  by the IPW method.
  In the comparison for time to event models, it is possible to use
  directly when \code{UP} and \code{DOWN} subjects are known.
}
\usage{
nricens.ipw.main(time, event, upp, dwn, t0)
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
}
