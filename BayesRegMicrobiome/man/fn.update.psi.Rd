\name{fn.update.psi}
\alias{fn.update.psi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	fn.update.psi
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.update.psi(Delta_1, DD, LL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Delta_1}{
%%     ~~Describe \code{Delta_1} here~~
}
  \item{DD}{
%%     ~~Describe \code{DD} here~~
}
  \item{LL}{
%%     ~~Describe \code{LL} here~~
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
function (Delta_1, DD, LL) 
{
    L_cnt <- table(c(Delta_1, (1:LL))) - 1
    psi <- rgamma(LL, L_cnt + DD, 1)
    psi <- psi/sum(psi)
    return(psi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
