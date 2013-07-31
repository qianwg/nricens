\name{nribin.count.main}
\alias{nribin.count.main}
\title{
  Estimate NRI by the counting method.
}
\description{
  Internaly used function by \code{nribin}.
  In the comparison for binary models, it is possible to use directly
  when \code{UP} and \code{DOWN} subjects are known.
}
\usage{
nribin.count.main(event, upp, dwn)
}
\arguments{
  \item{event}{
    Vector of event indicators, 1 for event of interest, 0 for non-event.
  }
  \item{upp, dwn}{
    Vector of logical values to determine subjects who belong to
    \code{UP} and \code{DOWN}, respectively.
  }
}
