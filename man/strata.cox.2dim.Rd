\name{strata.cox.2dim}
\alias{strata.cox.2dim}

\title{
  Determining subgroup based on two scores from the models.
}
\description{
  Determining subgroup based on a score from the working Cox models,
  for event and censoring time. Principal component analysis will be
  conducted for constructing subgoups.
  When one of the two working models are failed to obtain,
  the function \code{strata.cox.1dim} will be used for detemining subgroup.
  The total number of subgoups should be given.
}
\usage{
strata.cox.2dim(time, event, zcT, zcC, ms.strt, pca2nd, subs = NULL)
}

\arguments{
  \item{time}{
    Vector of time to event
}
  \item{event}{
    Vector of event indicator, 1 for event of interest, 0 for censoring
}
  \item{zcT}{
    Matrix of covariates for a working time-to-event model. Missing
    values are not allowed.
}
  \item{zcC}{
    Matrix of covariates for a working censoring time model. Missing
    values are not allowed.
}
  \item{ms.strt}{
    Scalar to define the number of minimum samples in a strata.
}
  \item{pca2nd}{
    Scalar for the number of strata of the second principle component.
}
  \item{subs}{
    Vector of Logical value to determine which subjects to use for a
    calculation. Default value is \code{NULL}, in which all subjects
    are used.
}
}
