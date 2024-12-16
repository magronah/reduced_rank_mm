\name{fn.initialize}
\alias{fn.initialize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	fn.initialize
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.initialize(hpara, YY, ZZ_1, Tot_N, PP, JJ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpara}{
%%     ~~Describe \code{hpara} here~~
}
  \item{YY}{
%%     ~~Describe \code{YY} here~~
}
  \item{ZZ_1}{
%%     ~~Describe \code{ZZ_1} here~~
}
  \item{Tot_N}{
%%     ~~Describe \code{Tot_N} here~~
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
function (hpara, YY, ZZ_1, Tot_N, PP, JJ) 
{
    cur_sam <- NULL
    sum_Y <- apply(YY, 1, sum)
    cur_sam$r <- sum_Y/sum(sum_Y)
    cur_sam$s <- rgamma(JJ, hpara$a_s, hpara$b_s)
    cur_sam$delta <- array(NA, dim = c(Tot_N, 2))
    cur_sam$delta[, 1] <- sample((1:hpara$L), Tot_N, TRUE, prob = rep(1, 
        hpara$L))
    cur_sam$delta[, 2] <- sample((0:1), Tot_N, TRUE, prob = c(1, 
        1))
    cur_sam$psi <- rep(1, hpara$L)/hpara$L
    cur_sam$w <- rbeta(hpara$L, hpara$a_w, hpara$b_w)
    cur_sam$eta <- rnorm(hpara$L, hpara$c_r, sqrt(hpara$v2_eta))
    cur_sam$beta <- array(0.01, dim = c(PP, JJ))
    cur_sam$beZ <- ZZ_1[, -(1:2)] \%*\% cur_sam$beta
    cur_sam$sig2 <- 1/rgamma(JJ, hpara$a_sig, hpara$b_sig)
    cur_sam$lam2 <- rgamma(JJ, hpara$a_lam, hpara$b_lam)
    cur_sam$phi <- array(NA, dim = c(PP, JJ))
    for (j in 1:JJ) {
        cur_sam$phi[, j] <- rgamma(PP, 1, cur_sam$lam2[j]/2)
    }
    r_mat <- rep.col(cur_sam$r, JJ)
    cur_sam$th0 <- log(apply((YY/r_mat), 2, median) + 0.01)
    cur_sam$delta_th0 <- array(NA, dim = c(JJ, 2))
    cur_sam$delta_th0[, 1] <- sample((1:hpara$L_th0), JJ, TRUE, 
        prob = rep(1, hpara$L_th0))
    cur_sam$delta_th0[, 2] <- sample((0:1), JJ, TRUE, prob = c(1, 
        1))
    cur_sam$psi_th0 <- rep(1, hpara$L_th0)/hpara$L_th0
    cur_sam$w_th0 <- rbeta(hpara$L_th0, hpara$aw_th0, hpara$bw_th0)
    cur_sam$eta_th0 <- rnorm(hpara$L_th0, hpara$c_th0, sqrt(hpara$v2_th0))
    cur_sam$th <- array(0, dim = c(hpara$M, JJ))
    cur_sam$Kth <- array(0, dim = c(Tot_N, JJ))
    th0_mat <- rep.row(cur_sam$th0, Tot_N)
    cur_sam$mu <- exp(th0_mat + cur_sam$Kth + cur_sam$beZ)
    cur_sam$tau2 <- 1/rgamma(JJ, hpara$a_tau, hpara$b_tau)
    return(cur_sam)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
