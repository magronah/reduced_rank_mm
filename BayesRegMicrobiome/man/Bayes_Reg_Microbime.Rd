\name{Bayes_Reg_Microbime}
\alias{Bayes_Reg_Microbime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayes_Reg_Microbime
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Bayes_Reg_Microbime(hpara, Dat, n_burn, n_sam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpara}{
%%     ~~Describe \code{hpara} here~~
}
  \item{Dat}{
%%     ~~Describe \code{Dat} here~~
}
  \item{n_burn}{
%%     ~~Describe \code{n_burn} here~~
}
  \item{n_sam}{
%%     ~~Describe \code{n_sam} here~~
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
function (hpara, Dat, n_burn, n_sam) 
{
    ini_sam <- fn.initialize(hpara, Dat$Y, Dat$Z_1, Dat$Tot_N, 
        Dat$P, Dat$J)
    print(date())
    set.seed(79861)
    ini_sam_1 <- fn.run.MCMC(hpara, Dat, ini_sam, n_burn)
    set.seed(79861)
    print(date())
    Save_MCMC_sam <- fn.run.MCMC.1(hpara, Dat, ini_sam_1, n_sam)
    return(Save_MCMC_sam)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
