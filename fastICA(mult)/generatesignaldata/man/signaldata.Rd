\name{signaldata}
\alias{signaldata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
signaldata(func, p, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
%%     ~~Describe \code{func} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (func, p, n) 
{
    S <- matrix(, p, 1)
    i = 0
    while (i < n) {
        if (func == "unif") {
            u <- runif(p, -sqrt(3), sqrt(3))
        }
        else if (func == "laplace") {
            u <- rlaplace(p, location = 0, scale = 1)
        }
        S <- cbind(S, u)
        i <- i + 1
    }
    S <- S[, -1]
    return(S)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
