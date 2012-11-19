\name{nricens.ipw.main}
\alias{nricens.ipw.main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nricens.ipw.main(time, event, upp, dwn, t0, strt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
%%     ~~Describe \code{time} here~~
}
  \item{event}{
%%     ~~Describe \code{event} here~~
}
  \item{upp}{
%%     ~~Describe \code{upp} here~~
}
  \item{dwn}{
%%     ~~Describe \code{dwn} here~~
}
  \item{t0}{
%%     ~~Describe \code{t0} here~~
}
  \item{strt}{
%%     ~~Describe \code{strt} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (time, event, upp, dwn, t0, strt) 
{
    case = time <= t0 & event == 1
    ctrl = time > t0
    w.case = rep(NA, length(time))
    w.ctrl = rep(NA, length(time))
    unq = unique(strt)
    for (i in unq) {
        ff = strt == i
        wk = get.surv.km(time[ff], as.i(event[ff] == 0), c(time[ff], 
            t0))
        w.case[ff] = 1/wk[-length(wk)]
        w.ctrl[ff] = 1/wk[length(wk)]
    }
    w.case[!case] = 0
    w.ctrl[!ctrl] = 0
    pr.upp.case = sum(w.case * upp)/sum(w.case)
    pr.dwn.case = sum(w.case * dwn)/sum(w.case)
    pr.dwn.ctrl = sum(w.ctrl * dwn)/sum(w.ctrl)
    pr.upp.ctrl = sum(w.ctrl * upp)/sum(w.ctrl)
    nri.case = pr.upp.case - pr.dwn.case
    nri.ctrl = pr.dwn.ctrl - pr.upp.ctrl
    nri = nri.case + nri.ctrl
    return(c(nri, nri.case, nri.ctrl, pr.upp.case, pr.dwn.case, 
        pr.dwn.ctrl, pr.upp.ctrl))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
