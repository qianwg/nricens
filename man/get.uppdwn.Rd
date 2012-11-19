\name{get.uppdwn}
\alias{get.uppdwn}

\title{
Obtaining UP and DOWN in the NRI calculation
}
\description{
This function calculates the UP and DOWN which is required for the NRI
calculation. The risk are calculated by the Cox regression model.
Two types of UP and DOWN definitions can be calculated: risk cutoff and
risk difference.
}
\usage{
get.uppdwn(time, event, z.std, z.new, t0, cut, region.method, msg = F)
}
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
    time of defining evnet/non-event.
}
  \item{cut}{
    Scalar or vector value.
    For a risk-cutoff NRI calculation, \code{cut} is a scalar or vector
    value to define risk cutoff value. For a risk-difference NRI
    calculation, \code{cut} is a scalar to define the threshold of the
    difference in risk. When \code{cut=0} is specified, continuous NRI
    will be provided. See also \code{region.method}.
}
  \item{region.method}{
    Character to define the region for detemining \code{UP} and
    \code{DOWN} in NRI calculation. For \code{"region.method = 'cutoff'"},
    the risk cutoff based NRI is calculated. For \code{"region.method = 'diff'"},
    the risk difference based NRI is calculated.
}
  \item{msg}{
    Logical value. If \code{msg = TRUE}, some information will be
    provided. Default value is \code{FALSE}.
}
}
