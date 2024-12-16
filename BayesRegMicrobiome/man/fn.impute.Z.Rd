\name{fn.impute.Z}
\alias{fn.impute.Z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	fn.impute.Z
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.impute.Z(YY, NN, Tot_N, JJ, Rep_n, Miss_ind, miss_cov_ind, Beta, Th0, KTh, RR, SS, Max_Cat, ZZ_1_cur, ZZ_cur, BeZ_cur, MMu_cur)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{YY}{
%%     ~~Describe \code{YY} here~~
}
  \item{NN}{
%%     ~~Describe \code{NN} here~~
}
  \item{Tot_N}{
%%     ~~Describe \code{Tot_N} here~~
}
  \item{JJ}{
%%     ~~Describe \code{JJ} here~~
}
  \item{Rep_n}{
%%     ~~Describe \code{Rep_n} here~~
}
  \item{Miss_ind}{
%%     ~~Describe \code{Miss_ind} here~~
}
  \item{miss_cov_ind}{
%%     ~~Describe \code{miss_cov_ind} here~~
}
  \item{Beta}{
%%     ~~Describe \code{Beta} here~~
}
  \item{Th0}{
%%     ~~Describe \code{Th0} here~~
}
  \item{KTh}{
%%     ~~Describe \code{KTh} here~~
}
  \item{RR}{
%%     ~~Describe \code{RR} here~~
}
  \item{SS}{
%%     ~~Describe \code{SS} here~~
}
  \item{Max_Cat}{
%%     ~~Describe \code{Max_Cat} here~~
}
  \item{ZZ_1_cur}{
%%     ~~Describe \code{ZZ_1_cur} here~~
}
  \item{ZZ_cur}{
%%     ~~Describe \code{ZZ_cur} here~~
}
  \item{BeZ_cur}{
%%     ~~Describe \code{BeZ_cur} here~~
}
  \item{MMu_cur}{
%%     ~~Describe \code{MMu_cur} here~~
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
function (YY, NN, Tot_N, JJ, Rep_n, Miss_ind, miss_cov_ind, Beta, 
    Th0, KTh, RR, SS, Max_Cat, ZZ_1_cur, ZZ_cur, BeZ_cur, MMu_cur) 
{
    Th0_Kth <- rep.row(Th0, Tot_N) + KTh
    RR_mat <- rep.col(RR, JJ)
    Z_impute <- array(-1, dim = c(NN, length(miss_cov_ind)))
    for (i_p in miss_cov_ind) {
        miss_t_i <- (1:NN)[Miss_ind[, (i_p + 2)] == 1]
        if (i_p > 1) {
            tmp <- sum(Max_Cat[1:(i_p - 1)])
        }
        else {
            tmp <- 0
        }
        cov_ind <- ((tmp + 1):(tmp + Max_Cat[i_p]))
        ZZ_cur[miss_t_i, 2 + cov_ind] <- 0
        for (t_i in miss_t_i) {
            ind_sam <- (1:Tot_N)[ZZ_1_cur[, 1] == t_i]
            ZZ_1_cur[ind_sam, 2 + cov_ind] <- 0
            Y_i <- YY[ind_sam, ]
            SS_mat <- rep.row(SS, Rep_n[t_i])
            BeZ_mat <- array(NA, dim = c(Rep_n[t_i], JJ, Max_Cat[i_p] + 
                1))
            mu_mat <- array(NA, dim = c(Rep_n[t_i], JJ, Max_Cat[i_p] + 
                1))
            ZZ_i <- as.matrix(ZZ_cur[t_i, -(1:2)])
            Prob_Z <- rep(NA, Max_Cat[i_p] + 1)
            i_c <- 0
            ZZ_i_tmp <- ZZ_i
            BeZ_mat[, , i_c + 1] <- rep.row((ZZ_i_tmp) \%*\% Beta, 
                Rep_n[t_i])
            mu_mat[, , i_c + 1] <- mu_tmp <- exp(Th0_Kth[ind_sam, 
                ] + BeZ_mat[, , i_c + 1])
            m_tmp <- RR_mat[ind_sam, ] * SS_mat * mu_tmp
            Prob_Z[i_c + 1] <- sum(Y_i * log(m_tmp/(1 + m_tmp)) - 
                (1/SS_mat) * log(1 + m_tmp))
            for (i_c in 1:Max_Cat[i_p]) {
                ZZ_i_tmp <- ZZ_i
                ZZ_i_tmp[cov_ind[i_c]] <- 1
                BeZ_mat[, , i_c + 1] <- rep.row((ZZ_i_tmp) \%*\% 
                  Beta, Rep_n[t_i])
                mu_mat[, , i_c + 1] <- mu_tmp <- exp(Th0_Kth[ind_sam, 
                  ] + BeZ_mat[, , i_c + 1])
                m_tmp <- RR_mat[ind_sam, ] * SS_mat * mu_tmp
                Prob_Z[i_c + 1] <- sum(Y_i * log(m_tmp/(1 + m_tmp)) - 
                  (1/SS_mat) * log(1 + m_tmp))
            }
            Prob_Z <- exp(Prob_Z - max(Prob_Z))
            c_sam <- sample((0:Max_Cat[i_p]), 1, FALSE, Prob_Z)
            if (c_sam > 0) {
                ZZ_cur[t_i, 2 + cov_ind[c_sam]] <- 1
                ZZ_1_cur[ind_sam, 2 + cov_ind[c_sam]] <- 1
            }
            Z_impute[t_i, i_p] <- c_sam
            MMu_cur[ind_sam, ] <- mu_mat[, , c_sam + 1]
            BeZ_cur[ind_sam, ] <- BeZ_mat[, , c_sam + 1]
        }
    }
    return(list(Z_1 = ZZ_1_cur, Z = ZZ_cur, BeZ = BeZ_cur, MMu = MMu_cur, 
        Impute_Z = Z_impute))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
