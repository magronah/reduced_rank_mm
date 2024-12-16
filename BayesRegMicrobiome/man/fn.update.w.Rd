\name{fn.update.w}
\alias{fn.update.w}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	fn.update.w
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.update.w(Delta, LL, aa_w, bb_w, Eta, RRt, C_R, u2, w_cur)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Delta}{
%%     ~~Describe \code{Delta} here~~
}
  \item{LL}{
%%     ~~Describe \code{LL} here~~
}
  \item{aa_w}{
%%     ~~Describe \code{aa_w} here~~
}
  \item{bb_w}{
%%     ~~Describe \code{bb_w} here~~
}
  \item{Eta}{
%%     ~~Describe \code{Eta} here~~
}
  \item{RRt}{
%%     ~~Describe \code{RRt} here~~
}
  \item{C_R}{
%%     ~~Describe \code{C_R} here~~
}
  \item{u2}{
%%     ~~Describe \code{u2} here~~
}
  \item{w_cur}{
%%     ~~Describe \code{w_cur} here~~
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
function (Delta, LL, aa_w, bb_w, Eta, RRt, C_R, u2, w_cur) 
{
    Delta1 <- Delta[, 1]
    Delta2 <- Delta[, 2]
    tmp <- log(w_cur/(1 - w_cur))
    tmp <- tmp + rnorm(LL, 0, 0.5)
    w_pro <- exp(tmp)/(1 + exp(tmp))
    for (l in 1:LL) {
        w_l_cur <- w_cur[l]
        w_l_pro <- w_pro[l]
        Delta_l0 <- (Delta1 == l) * (Delta2 == 0)
        Delta_l1 <- (Delta1 == l) * (Delta2 == 1)
        RRt_l <- RRt[Delta_l1]
        aa_1 <- aa_w + sum(Delta_l0)
        bb_1 <- bb_w + sum(Delta_l1)
        m_cur <- (C_R - w_l_cur * Eta[l])/(1 - w_l_cur)
        m_pro <- (C_R - w_l_pro * Eta[l])/(1 - w_l_pro)
        A_cur <- aa_1 * log(w_l_cur) + (bb_1) * log(1 - w_l_cur) - 
            sum((RRt_l - m_cur)^2)/2/u2
        A_pro <- aa_1 * log(w_l_pro) + (bb_1) * log(1 - w_l_pro) - 
            sum((RRt_l - m_pro)^2)/2/u2
        if (log(runif(1)) < (A_pro - A_cur)) {
            w_cur[l] <- w_l_pro
        }
    }
    return(w_cur)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
