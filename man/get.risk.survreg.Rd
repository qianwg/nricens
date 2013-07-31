\name{get.risk.survreg}
\alias{get.risk.survreg}
\title{
  Calculate individual risks based on the parametric survival regression
  model.
}
\description{
  This function estimates \eqn{Pr\{T < t_0\}} given covariates under the
  parametric survival model.
}
\usage{
get.risk.survreg(mdl, t0)
}
\arguments{
  \item{mdl}{
    \code{survreg} object.
  }
  \item{t0}{
    Scalar value indicating a time to determine evnet/non-event (\eqn{t_0}).
  }
}
