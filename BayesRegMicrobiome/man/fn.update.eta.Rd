\name{fn.update.eta}
\alias{fn.update.eta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	fn.update.eta
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.update.eta(Delta, LL, ww, u2, v2_eta, RRt, C_R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Delta}{
%%     ~~Describe \code{Delta} here~~
}
  \item{LL}{
%%     ~~Describe \code{LL} here~~
}
  \item{ww}{
%%     ~~Describe \code{ww} here~~
}
  \item{u2}{
%%     ~~Describe \code{u2} here~~
}
  \item{v2_eta}{
%%     ~~Describe \code{v2_eta} here~~
}
  \item{RRt}{
%%     ~~Describe \code{RRt} here~~
}
  \item{C_R}{
%%     ~~Describe \code{C_R} here~~
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
function (Delta, LL, ww, u2, v2_eta, RRt, C_R) 
{
    w_1_w <- ww/(1 - ww)
    Cr_1_w <- C_R/(1 - ww)
    Delta1 <- Delta[, 1]
    Delta2 <- Delta[, 2]
    Eta <- rep(NA, LL)
    for (l in 1:LL) {
        ind_l0 <- ((Delta1 == l) * (Delta2 == 0))
        ind_l1 <- ((Delta1 == l) * (Delta2 == 1))
        v2_tmp <- 1/v2_eta + sum(ind_l0)/u2 + sum(ind_l1) * (w_1_w[l])^2/u2
        v2_tmp <- 1/v2_tmp
        m_tmp <- C_R/v2_eta + (sum(ind_l0 * RRt) - w_1_w[l] * 
            sum(ind_l1 * (RRt - Cr_1_w[l])))/u2
        m_tmp <- v2_tmp * m_tmp
        Eta[l] <- rnorm(1, m_tmp, sqrt(v2_tmp))
    }
    return(Eta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
