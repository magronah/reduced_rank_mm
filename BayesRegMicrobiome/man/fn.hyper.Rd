\name{fn.hyper}
\alias{fn.hyper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	fn.hyper
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.hyper(YY, Tot_N, std_t, TT, PP, JJ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{YY}{
%%     ~~Describe \code{YY} here~~
}
  \item{Tot_N}{
%%     ~~Describe \code{Tot_N} here~~
}
  \item{std_t}{
%%     ~~Describe \code{std_t} here~~
}
  \item{TT}{
%%     ~~Describe \code{TT} here~~
}
  \item{PP}{
%%     ~~Describe \code{PP} here~~
}
  \item{JJ}{
%%     ~~Describe \code{JJ} here~~
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
function (YY, Tot_N, std_t, TT, PP, JJ) 
{
    hpara <- NULL
    hpara$a_s <- 1
    hpara$b_s <- 2
    hpara$L <- 30
    hpara$u2 <- 1
    sum_Y <- apply(YY, 1, sum)
    r_mle <- (sum_Y/sum(sum_Y))
    hpara$c_r <- mean(log(r_mle))
    hpara$a_w <- 1
    hpara$b_w <- 1
    hpara$v2_eta <- 1
    hpara$d <- 10
    hpara$L_th0 <- 50
    hpara$u2_th0 <- 1
    r_mle <- rep.col(r_mle, JJ)
    th0_mle <- log(apply(YY/r_mle, 2, mean))
    hpara$c_th0 <- mean(th0_mle)
    hpara$aw_th0 <- 1
    hpara$bw_th0 <- 1
    hpara$v2_th0 <- 2
    hpara$d_th0 <- 10
    hpara$a_sig <- 3
    hpara$b_sig <- 3
    hpara$a_lam <- 0.5
    hpara$b_lam <- 0.5
    hpara$M <- 13
    tt0 <- 10
    hpara$u <- seq(-tt0, TT + tt0, length.out = hpara$M)
    hpara$gam2 <- ((TT + tt0 * 2)/hpara$M)^2
    hpara$K <- array(NA, dim = c(Tot_N, hpara$M))
    for (i_m in 1:hpara$M) {
        hpara$K[, i_m] <- dnorm(std_t, hpara$u[i_m], sqrt(hpara$gam2))
    }
    hpara$a_tau <- 0.1
    hpara$b_tau <- 0.1
    return(hpara)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
