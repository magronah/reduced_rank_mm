\name{fn.run.MCMC.1}
\alias{fn.run.MCMC.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	fn.run.MCMC.1
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn.run.MCMC.1(hpara, Dat, cur_sam, N_iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpara}{
%%     ~~Describe \code{hpara} here~~
}
  \item{Dat}{
%%     ~~Describe \code{Dat} here~~
}
  \item{cur_sam}{
%%     ~~Describe \code{cur_sam} here~~
}
  \item{N_iter}{
%%     ~~Describe \code{N_iter} here~~
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
function (hpara, Dat, cur_sam, N_iter) 
{
    YY <- Dat$Y
    JJ <- Dat$J
    PP <- Dat$P
    Tot_N <- Dat$Tot_N
    NN <- Dat$n
    Rep_n <- Dat$rep_K
    Miss_ind <- Dat$Miss_ind
    miss_cov_ind <- Dat$miss_cov_ind
    miss_sam_ind <- (1:NN)[apply(Miss_ind[, -(1:2)] == 1, 1, 
        sum) > 0]
    N_miss_sam <- length(miss_sam_ind)
    Impute_Z <- array(-1, dim = c(N_miss_sam, length(miss_cov_ind)))
    MM <- hpara$M
    YY_c <- array(YY, dim = c(Tot_N * JJ, 1))[, 1]
    KK_c <- array(hpara$K, dim = c(Tot_N * hpara$M, 1))[, 1]
    ZZ <- cur_sam$Z
    ZZ_1 <- cur_sam$Z_1
    ZZ_1_cov <- cur_sam$Z_1_cov
    cur_sam$Z <- NULL
    cur_sam$Z_1 <- NULL
    cur_sam$Z_1_cov <- NULL
    SAVE_sam <- NULL
    SAVE_sam$beta <- array(NA, dim = c(PP, JJ, N_iter))
    SAVE_sam$phi <- array(NA, dim = c(PP, JJ, N_iter))
    SAVE_sam$s <- array(NA, dim = c(JJ, N_iter))
    SAVE_sam$r <- array(NA, dim = c(Tot_N, N_iter))
    SAVE_sam$th0 <- array(NA, dim = c(JJ, N_iter))
    SAVE_sam$th <- array(NA, dim = c(hpara$M, JJ, N_iter))
    SAVE_sam$lam2 <- SAVE_sam$tau2 <- SAVE_sam$sig2 <- array(NA, 
        dim = c(JJ, N_iter))
    SAVE_sam$Z_1_cov <- array(NA, dim = c(Tot_N, 10, N_iter))
    print("after burn-in")
    for (i_iter in 1:N_iter) {
        if ((i_iter\%\%1000) == 0) {
            print(paste("i.iter=", i_iter))
            print(date())
        }
        cur_sam$delta <- fn.update.delta(cur_sam$delta, log(cur_sam$psi), 
            cur_sam$w, cur_sam$eta, hpara$c_r, log(cur_sam$r), 
            hpara$u2, Tot_N, hpara$L)
        cur_sam$eta <- fn.update.eta(cur_sam$delta, hpara$L, 
            cur_sam$w, hpara$u2, hpara$v2_eta, log(cur_sam$r), 
            hpara$c_r)
        cur_sam$w <- fn.update.w(cur_sam$delta, hpara$L, hpara$a_w, 
            hpara$b_w, cur_sam$eta, log(cur_sam$r), hpara$c_r, 
            hpara$u2, cur_sam$w)
        cur_sam$psi <- fn.update.psi(cur_sam$delta[, 1], hpara$d, 
            hpara$L)
        cur_sam$delta_th0 <- fn.update.delta(cur_sam$delta_th0, 
            log(cur_sam$psi_th0), cur_sam$w_th0, cur_sam$eta_th0, 
            hpara$c_th0, cur_sam$th0, hpara$u2_th0, JJ, hpara$L_th0)
        cur_sam$eta_th0 <- fn.update.eta(cur_sam$delta_th0, hpara$L_th0, 
            cur_sam$w_th0, hpara$u2_th0, hpara$v2_th0, cur_sam$th0, 
            hpara$c_th0)
        cur_sam$w_th0 <- fn.update.w(cur_sam$delta_th0, hpara$L_th0, 
            hpara$aw_th0, hpara$bw_th0, cur_sam$eta_th0, cur_sam$th0, 
            hpara$c_th0, hpara$u2_th0, cur_sam$w_th0)
        cur_sam$psi_th0 <- fn.update.psi(cur_sam$delta_th0[, 
            1], hpara$d_th0, hpara$L_th0)
        Th_c <- array(cur_sam$th, dim = c(MM * JJ, 1))[, 1]
        KTh_c <- array(cur_sam$Kth, dim = c(Tot_N * JJ, 1))[, 
            1]
        Beta_c <- array(cur_sam$beta, dim = c(JJ * PP, 1))[, 
            1]
        BeZ_c <- array(cur_sam$beZ, dim = c(Tot_N * JJ, 1))[, 
            1]
        MMu_c <- array(cur_sam$mu, dim = c(Tot_N * JJ, 1))[, 
            1]
        Phi_c <- array(cur_sam$phi, dim = c(JJ * PP, 1))[, 1]
        ZZ_c <- array(ZZ_1_cov, dim = c(Tot_N * PP, 1))[, 1]
        Rt_mean <- c(cur_sam$eta, (hpara$c_r - cur_sam$w * cur_sam$eta)/(1 - 
            cur_sam$w))
        ind <- cur_sam$delta[, 1] + hpara$L * cur_sam$delta[, 
            2]
        Rt_mean <- Rt_mean[ind]
        Th0_mean <- c(cur_sam$eta_th0, (hpara$c_th0 - cur_sam$w_th0 * 
            cur_sam$eta_th0)/(1 - cur_sam$w_th0))
        ind <- cur_sam$delta_th0[, 1] + hpara$L_th0 * cur_sam$delta_th0[, 
            2]
        Th0_mean <- Th0_mean[ind]
        output <- .C("fn_update_Beta_Theta", JJ = as.integer(JJ), 
            PP = as.integer(PP), MM = as.integer(MM), Tot_N = as.integer(Tot_N), 
            YY = as.double(YY_c), a_sig = as.double(hpara$a_sig), 
            b_sig = as.double(hpara$b_sig), sig2 = as.double(cur_sam$sig2), 
            a_tau = as.double(hpara$a_tau), b_tau = as.double(hpara$b_tau), 
            tau2 = as.double(cur_sam$tau2), a_lam = as.double(hpara$a_lam), 
            b_lam = as.double(hpara$b_lam), lam2 = as.double(cur_sam$lam2), 
            Rt_m = as.double(Rt_mean), u2 = as.double(hpara$u2), 
            RR = as.double(cur_sam$r), a_s = as.double(hpara$a_s), 
            b_s = as.double(hpara$b_s), SS = as.double(cur_sam$s), 
            Th0_bar = as.double(Th0_mean), u2_th0 = as.double(hpara$u2_th0), 
            Th0 = as.double(cur_sam$th0), KK = as.double(KK_c), 
            Th = as.double(Th_c), KTh = as.double(KTh_c), ZZ = as.double(ZZ_c), 
            phi = as.double(Phi_c), Beta = as.double(Beta_c), 
            BeZ = as.double(BeZ_c), MMu = as.double(MMu_c))
        cur_sam$th0 <- output$Th0
        cur_sam$th <- array(output$Th, dim = c(MM, JJ))
        cur_sam$Kth <- array(output$KTh, dim = c(Tot_N, JJ))
        cur_sam$beta <- array(output$Beta, dim = c(PP, JJ))
        cur_sam$beZ <- array(output$BeZ, dim = c(Tot_N, JJ))
        cur_sam$mu <- array(output$MMu, dim = c(Tot_N, JJ))
        cur_sam$s <- output$SS
        cur_sam$tau2 <- output$tau2
        cur_sam$sig2 <- output$sig2
        cur_sam$lam2 <- output$lam2
        cur_sam$phi <- array(output$phi, dim = c(PP, JJ))
        cur_sam$r <- output$RR
        tmp <- fn.impute.Z(YY, NN, Tot_N, JJ, Rep_n, Miss_ind, 
            miss_cov_ind, cur_sam$beta, cur_sam$th0, cur_sam$Kth, 
            cur_sam$r, cur_sam$s, Dat$max_cat, ZZ_1, ZZ, cur_sam$beZ, 
            cur_sam$mu)
        cur_sam$mu <- tmp$MMu
        cur_sam$beZ <- tmp$BeZ
        ZZ <- tmp$Z
        ZZ_1 <- tmp$Z_1
        ZZ_1_cov <- ZZ_1[, -c(1:2)]
        Impute_Z <- tmp$Impute_Z
        SAVE_sam$beta[, , i_iter] <- cur_sam$beta
        SAVE_sam$phi[, , i_iter] <- cur_sam$phi
        SAVE_sam$s[, i_iter] <- cur_sam$s
        SAVE_sam$r[, i_iter] <- cur_sam$r
        SAVE_sam$th0[, i_iter] <- cur_sam$th0
        SAVE_sam$th[, , i_iter] <- cur_sam$th
        SAVE_sam$lam2[, i_iter] <- cur_sam$lam2
        SAVE_sam$sig2[, i_iter] <- cur_sam$sig2
        SAVE_sam$tau2[, i_iter] <- cur_sam$tau2
        SAVE_sam$Z_1_cov[, , i_iter] <- ZZ_1_cov[, 1:10]
    }
    return(SAVE_sam)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
