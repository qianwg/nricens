\name{get.uppdwn.bin}
\alias{get.uppdwn.bin}
\title{
  Determine \code{UP} and \code{DOWN} for the NRI calculation
}
\description{
  Internaly used function to detemine subjects who belong to \code{UP}
  and \code{DOWN}.
}
\usage{
get.uppdwn.bin(event, objs, flag.mdl, flag.prd, flag.rsk, updown, cut,
               link, msg = FALSE)
}

\arguments{
  \item{event}{
    Vector of event indicators, 1 for event of interest, 0 for non-event.
  }
  \item{objs}{
    List of data.
  }
  \item{flag.mdl, flag.prd, flag.rsk}{
    Logical values to determine the type of calculation based on the
    input data.
  }
  \item{updown}{
    Character to specify the method to determine \code{UP} and
    \code{DOWN}.
  }
  \item{cut}{
    Scalar or vector to specify the cutoff value(s) of predicted risks
    for determining \code{UP} and \code{DOWN}.
    See also \code{cut} in \code{nricens} or \code{nribin}.
  }
  \item{link}{
    Character to specify the link function for a glm fitting.
  }
  \item{msg}{
    Logical value to display computation process.
    Setting \code{FALSE} leads silent execution.
  }
}
