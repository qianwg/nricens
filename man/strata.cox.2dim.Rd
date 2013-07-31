\name{strata.cox.2dim}
\alias{strata.cox.2dim}

\title{
  Determine strata for the weighted Kaplan-Meier estimator.
}
\description{
  Determine strata for estimating a survival function by the weighted
  Kaplan-Meier estimator (Hsu and Taylor, 2010).
  Strata are constructed based on scores from working Cox models for an
  event and a censoring time.
  Stratification is conducted using the principal components of these two
  scores.
}
\usage{
strata.cox.2dim(time, event, zc.T, zc.C, m, pca2nd, subs = NULL)
}
\arguments{
  \item{time}{
    Vector of follow up times.
  }
  \item{event}{
    Vector of event indicators, 1 for event of interest, 0 for censoring.
  }
  \item{zc.T, zc.C}{
    Matrix of covariates for working event and censoring time Cox
    models, respectively. Missing values are not allowed.
  }
  \item{m}{
    Scalar value to determine the number of minimum samples for each
    stratum.
  }
  \item{pca2nd}{
    Scalar value to determine the number of strata in the second
    principal component.
  }
  \item{subs}{
    Vector of logical values to determine which subjects to be used for
    a calculation. Default value is \code{NULL}, in this case all
    subjects are used for a calculation.
  }
}
\value{
  Returns a list of the following items:
  \item{strt}{
    Strata number for each subject.
  }
  \item{mdl.T, mdl.C}{
    Fitted \code{coxph} objects corresponding to an event and a censoring
    time model, respectively.
  }
  \item{pca}{
    Result of a principal component analysis by \code{prcomp}
    using scores of working Cox models.
  }
}
