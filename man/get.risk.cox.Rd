\name{get.risk.cox}
\alias{get.risk.cox}
\title{
  Calculating individual predicted risk by the Cox regression model.
}
\description{
  Calculating individual predicted risk by the Cox regression model.
}
\usage{
get.risk.cox(mdl, t0)
}
\arguments{
  \item{mdl}{
    Object derived from the \code{coxph}
  }
  \item{t0}{
    Time of defining evnet/non-event
  }
}
