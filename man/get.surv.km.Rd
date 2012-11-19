\name{get.surv.km}
\alias{get.surv.km}
\title{
Calculating the survival probability by the standard Kaplan-Meier
estimator at the given time.
}
\description{
Calculating the survival probability by the standard Kaplan-Meier
estimator given the time. \code{survfit} is run internally.
}
\usage{
get.surv.km(time, event, t0, subst = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
    Vector of time to event
}
  \item{event}{
    Vector of event indicator, 1 for event of interest, 0 for censoring
}
  \item{t0}{
    Time of defining evnet/non-event.
}
  \item{subst}{
    Vector of logical values to define subjects for a calculation.
    When \code{subst=NULL}, all subjects are used.
}
}
