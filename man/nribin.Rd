\name{nribin}
\alias{nribin}
\title{
NRI for binary models
}
\description{
  This function estimates the NRI for competing risk prediction models
  with binary response variable. \code{glm} objects, predictors, and
  predicted risks can be used as input data for the calculation.
  The risk category based NRI and the risk difference based NRI can be
  calculated.
  The percentile bootstrap method is used for an interval estimation.
} 
\usage{
nribin(event = NULL, mdl.std = NULL, mdl.new = NULL, z.std = NULL, z.new = NULL,
       p.std = NULL, p.new = NULL, updown = "category", cut = NULL,
       link = "logit", niter = 1000, alpha = 0.05, msg = TRUE)
       
}
\arguments{
  \item{event}{
    Vector of event indicators, 1 for event of interest, 0 for non-event.
  }
  \item{mdl.std, mdl.new}{
    \code{glm} objects corresponding to a standard and a new risk
    prediction models, respectively.
    Since predictors are extracted from these objects, \code{x=TRUE} is
    required for a glm fitting.
  }
  \item{z.std, z.new}{
    Matrix of predictors for a standard and a new risk prediction
    model, respectively.
    Neither factor nor character nor missing values are allowed.
  }
  \item{p.std, p.new}{
    Vector of predicted risks from a standard and a new prediction
    model, respectively.
  }
  \item{updown}{
    Character to specify the method to determine \code{UP} and
    \code{DOWN}.
    When \code{"category"} is specified (by default), the risk category
    based NRI is calculated.
    When \code{"diff"} is specified, the risk difference based NRI is
    calculated.
  }
  \item{cut}{
    Scalar or vector values to specify the cutoff value(s) of predicted
    risks for determining \code{UP} and \code{DOWN}.
    For the risk category based NRI (\code{updown = "category"}), this
    option corresponds to cutoff value(s) of risk categories.
    For the risk difference based NRI (\code{updown = "diff"}), this
    option corresponds to a cutoff value of a risk difference
    (only scalar is allowed).
  }
  \item{niter}{
    Scalar value to determine the number of bootstratp sampling. When
    \code{0} is specified, an interval estimation is skipped.
  }
  \item{link}{
    Character to specify a link function for a glm fitting.
  }
  \item{alpha}{
    1-alpha confidence intervals are calcualted.
  }
  \item{msg}{
    Logical value to display computation process.
    Setting \code{FALSE} leads silent execution.
  }
}
\details{
  Either one set of the following arguments should be specified for the
  NRI calculation: (\code{mdl.std}, \code{mdl.new});
  (\code{event}, \code{z.std}, \code{z.new});
  and (\code{event}, \code{p.std}, \code{p.new}).

  In the first set of the argument, (\code{mdl.std}, \code{mdl.new}),
  fitted results are used for the NRI calculation.
  \code{event}, \code{z.std}, and \code{z.new} are extracted from fitted
  result objects.
  The variance of model parameters are accounted for an interval
  estimation of the NRI.
  When \code{event} is specified in arguments, those specified is used
  without extracting from \code{glm} object.

  In the second set of the argument, (\code{event}, \code{z.std},
  \code{z.new}), a standard and a new prediction models are fitted
  inside this function with specified \code{link}.
  The variance of model parameters are also accounted for an interval
  estimation of the NRI.

  In the third set of the argument, (\code{event}, \code{p.std},
  \code{p.new}), predicted risks are used.
  Since fit of prediction models are not conducted while in a bootstrap,
  this can be used for a validation study by an external data source or
  by a cross-validation.

  For the risk category based NRI calculation, cutoff values of risk
  category can be specified by \code{cut}, which is a scalar for the
  case of two risk categories and is a vector for the case of more than
  two risk categories.
  \code{UP} and \code{DOWN} are determined by the movement in risk
  categories.

  For the risk difference based NRI calculation, cutoff values of risk
  difference can also be specified by \code{cut}, where \code{UP} and
  \code{DOWN} are defiend as \eqn{p_{new} - p_{standard} > \delta} and
  \eqn{p_{standard} - p_{new} > \delta}, respectively.
  \eqn{p_{standard}} and \eqn{p_{new}} are predicted individual risks
  from a standard and a new prediction model, respectively, and
  \eqn{\delta} corresponds to \code{cut}.
  The continuous NRI, which is the special version of the risk
  difference based NRI, can be calculated by specifying both
  \code{updown = "diff"} and \code{cut = 0}.

  Interval estimation is based on the percentile bootstrap method.
}

\value{
  Returns a list of the following items:
  \item{nri}{
    Point and interval estimates of the NRI and its components.
  }
  \item{mdl.std, mdl.new}{
    Fitted \code{glm} objects corresponding to a standard and a new
    prediction model, respectively.
    These items are provided when prediction models are fitted inside
    this function. Otherwise \code{NULL}.
  }
  \item{z.std, z.new}{
    Predictors of a standard and a new prediction model, respectively.
    These items are provided when they are extracted from \code{mdl.std}
    and \code{mdl.new}. Otherwise \code{NULL}.
  }
  \item{p.std, p.new}{
    Predicted risks by a standard and a new prediction model,
    respectively.
  }
  \item{up, down}{
    Logical values to show subjects who belong to \code{UP} and
    \code{DOWN}, respectively.
  }
  \item{rtab, rtab.case, rtab.ctrl}{
    \code{table} objects corresponding to reclassification tables for all,
    case, and control subjects, respectively.
    These items are provided only when the risk category based NRI is
    calculated and \code{msg = TRUE}.
  }
  \item{bootstrapsample}{
    Results of each bootstrap sample.
  }
}
\examples{
## here consider pbc dataset in survival package as an example
library(survival)
dat = pbc[1:312,]
dat$sex = ifelse(dat$sex=='f', 1, 0)

## subjects censored before 2000 days are excluded
dat = dat[ dat$time > 2000 | (dat$time < 2000 & dat$status == 2), ]

## predciting the event of 'death' before 2000 days
event = ifelse(dat$time < 2000 & dat$status == 2, 1, 0)

## standard prediction model: age, bilirubin, and albumin
z.std = as.matrix(subset(dat, select = c(age, bili, albumin)))

## new prediction model: age, bilirubin, albumin, and protime
z.new = as.matrix(subset(dat, select = c(age, bili, albumin, protime)))

## glm fit (logistic model)
mstd = glm(event ~ ., binomial(logit), data.frame(event, z.std), x=TRUE)
mnew = glm(event ~ ., binomial(logit), data.frame(event, z.new), x=TRUE)

## predicted risk
p.std = mstd$fitted.values
p.new = mnew$fitted.values

## Calculation of risk difference NRI using ('mdl.std', 'mdl.std').
nribin(mdl.std = mstd, mdl.new = mnew, cut = 0.02, niter = 0,
       updown = 'diff')

## Calculation of risk difference NRI using ('event', 'z.std', 'z.std').
nribin(event = event, z.std = z.std, z.new = z.new, cut = 0.02,
       niter = 0, updown = 'diff')

## Calculation of risk difference NRI using ('event', 'p.std', 'p.std').
nribin(event = event, p.std = p.std, p.new = p.new, cut = 0.02,
       niter = 0, updown = 'diff')

## Calculation of risk category NRI using ('mdl.std', 'mdl.std').
nribin(mdl.std = mstd, mdl.new = mnew, cut = c(0.2, 0.4),
       niter = 0, updown = 'category')
}
