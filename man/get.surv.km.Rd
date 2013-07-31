\name{get.surv.km}
\alias{get.surv.km}
\title{
  Calculate a survival probability at a given time by the standard
  Kaplan-Meier estimator.
}
\description{
  This function estimates \eqn{Pr\{T > t_0\}} by the Kaplan-Meier
  estimator. \eqn{t_0} should be given.
}
\usage{
get.surv.km(time, event, t0, subs = NULL)
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
  \item{subs}{
    Vector of logical values to determine which subjects to be used for
    a calculation. When this option is not specified
    (\code{subs = NULL}), all subjects are used for a calculation.
  }
}
