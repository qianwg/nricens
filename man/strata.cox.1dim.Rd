\name{strata.cox.1dim}
\alias{strata.cox.1dim}

\title{
  Determining subgroup based on a score from the model.
}
\description{
  Determining subgroup based on a score from the working Cox model,
  which is determined by the categorization of the linear predictor of
  working model. The number of subgoups should be given.
  When the score from the working model for the censoring time,
  censoring time should be treated as "event".
}
\usage{
strata.cox.1dim(time, event, zc, ms.strt)
}

\arguments{
  \item{time}{
    Vector of time to event
}
  \item{event}{
    Vector of event indicator, 1 for event of interest, 0 for censoring
}
  \item{zc}{
    Matrix of covariates for a working model.
    Missing values are not allowed.
}
  \item{ms.strt}{
    Scalar to define the number of minimum samples in a subgroup.
}
}
