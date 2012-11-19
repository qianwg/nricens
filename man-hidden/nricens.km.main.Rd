\name{nricens.km.main}
\alias{nricens.km.main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nricens.km.main(time, event, upp, dwn, t0, strt, strt.upp, strt.dwn)
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
  \item{strt.upp}{
%%     ~~Describe \code{strt.upp} here~~
}
  \item{strt.dwn}{
%%     ~~Describe \code{strt.dwn} here~~
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
function (time, event, upp, dwn, t0, strt, strt.upp, strt.dwn) 
{
    pr.upp = mean(upp)
    pr.dwn = mean(dwn)
    if (pr.upp == 0 && pr.dwn == 0) 
        return(rep(0, 7))
    if (!(length(unique(strt)) == 1 && length(unique(strt.upp)) == 
        1 && length(unique(strt.dwn)) == 1)) {
        if (pr.upp != 0) 
            strt[upp] = strt.upp[upp] + max(strt)
        if (pr.dwn != 0) 
            strt[dwn] = strt.dwn[dwn] + max(strt)
    }
    pr.ctrl = get.surv.wkm(time, as.i(event == 1), t0, strt)
    pr.ctrl.upp = pr.upp.ctrl = pr.upp.case = 0
    if (pr.upp != 0) {
        pr.ctrl.upp = get.surv.wkm(time[upp], as.i(event[upp] == 
            1), t0, strt.upp[upp])
        pr.upp.ctrl = pr.ctrl.upp * pr.upp/pr.ctrl
        pr.upp.case = (1 - pr.ctrl.upp) * pr.upp/(1 - pr.ctrl)
    }
    pr.ctrl.dwn = pr.dwn.ctrl = pr.dwn.case = 0
    if (pr.dwn != 0) {
        pr.ctrl.dwn = get.surv.wkm(time[dwn], as.i(event[dwn] == 
            1), t0, strt.dwn[dwn])
        pr.dwn.ctrl = pr.ctrl.dwn * pr.dwn/pr.ctrl
        pr.dwn.case = (1 - pr.ctrl.dwn) * pr.dwn/(1 - pr.ctrl)
    }
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
