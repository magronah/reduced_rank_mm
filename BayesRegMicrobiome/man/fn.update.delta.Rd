\name{fn.update.delta}
\alias{fn.update.delta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	fn.update.delta
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.update.delta(Delta_cur, log_Psi, ww, Eta, CR, RRt, u2, Tot_N, LL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Delta_cur}{
%%     ~~Describe \code{Delta_cur} here~~
}
  \item{log_Psi}{
%%     ~~Describe \code{log_Psi} here~~
}
  \item{ww}{
%%     ~~Describe \code{ww} here~~
}
  \item{Eta}{
%%     ~~Describe \code{Eta} here~~
}
  \item{CR}{
%%     ~~Describe \code{CR} here~~
}
  \item{RRt}{
%%     ~~Describe \code{RRt} here~~
}
  \item{u2}{
%%     ~~Describe \code{u2} here~~
}
  \item{Tot_N}{
%%     ~~Describe \code{Tot_N} here~~
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
function (Delta_cur, log_Psi, ww, Eta, CR, RRt, u2, Tot_N, LL) 
{
    P1 <- c((log_Psi + log(ww)), (log_Psi + log(1 - ww)))
    Mean_tmp <- c(Eta, (CR - ww * Eta)/(1 - ww))
    LL_grid <- (1:(2 * LL))
    for (i in 1:Tot_N) {
        Prob <- P1 - (RRt[i] - Mean_tmp)^2/2/u2
        Prob <- exp(Prob - max(Prob))
        ind <- sample(LL_grid, 1, FALSE, Prob)
        Delta2 <- 1 * (ind > LL)
        Delta1 <- ind - (Delta2 * LL)
        Delta_cur[i, ] <- c(Delta1, Delta2)
    }
    return(Delta_cur)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
