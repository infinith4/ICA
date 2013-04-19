\name{rlaplace}
\alias{rlaplace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rlaplace(n, location = 0, scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{location}{
%%     ~~Describe \code{location} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
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
function (n, location = 0, scale = 1) 
{
    if (!is.Numeric(n, positive = TRUE, integer.valued = TRUE, 
        allowable.length = 1)) 
        stop("bad input for argument 'n'")
    if (!is.Numeric(scale, positive = TRUE)) 
        stop("'scale' must be positive")
    location = rep(location, length.out = n)
    scale = rep(scale, length.out = n)
    r = runif(n)
    location - sign(r - 0.5) * scale * log(2 * ifelse(r < 0.5, 
        r, 1 - r))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
