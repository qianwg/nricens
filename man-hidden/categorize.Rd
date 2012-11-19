\name{categorize}
\alias{categorize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
categorize(dat, threshold = NULL, nlevel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{nlevel}{
%%     ~~Describe \code{nlevel} here~~
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
function (dat, threshold = NULL, nlevel = NULL) 
{
    if (is.null(threshold)) 
        threshold = quantile(dat, seq(0, 1, 1/nlevel), na.rm = T)[-1]
    else nlevel = length(threshold)
    threshold[nlevel] = max(dat, na.rm = T) + 1
    new = rep(NA, length(dat))
    for (i in 1:nlevel) new[is.na(new) & dat < threshold[i]] = i
    return(new)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
