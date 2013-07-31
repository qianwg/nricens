\name{get.risk.coxph}
\alias{get.risk.coxph}
\title{
  Calculate individual risks based on the Cox regression model.
}
\description{
  This function estimates \eqn{Pr\{T < t_0\}} given covariates under the
  Cox model.
}
\usage{
get.risk.coxph(mdl, t0)
}
\arguments{
  \item{mdl}{
    \code{coxph} object.
  }
  \item{t0}{
    Scalar value indicating a time to determine evnet/non-event (\eqn{t_0}).
  }
}
