\name{strata.1dim}
\alias{strata.1dim}
\title{
  Categorization of a continuous variable.
}
\description{
  Internaly used function to categorize a continuous variable.
}
\usage{
strata.1dim(score, m, subs = NULL)
}
\arguments{
  \item{score}{
    Numerical vector for a categorization.
  }
  \item{m}{
    Scalar value to determine the number of minimum samples for each
    stratum.
  }
  \item{subs}{
    Vector of logical values to determine which subjects to be used for
    a calculation. When this option is not specified (\code{"strt = NULL"}),
    all subjects are used for a calculation.
  }
}
