

Bayes_Reg_Microbime <- function(hpara, Dat, n_burn, n_sam)
{
    ## initialize random parameters
    ini_sam <- fn.initialize(hpara, Dat$Y, Dat$Z_1, Dat$Tot_N, Dat$P, Dat$J)
    
    ### run MCMC
    print(date())
    
    set.seed(79861)
    ini_sam_1 <- fn.run.MCMC(hpara, Dat, ini_sam, n_burn)
    
    ##
    set.seed(79861)
    print(date())
    Save_MCMC_sam <- fn.run.MCMC.1(hpara, Dat, ini_sam_1, n_sam)
    
    return(Save_MCMC_sam)
}





###  hpara <- hyper; Dat <- SIM_dat; cur_sam <- ini_sam; N_iter <- N_burn
### hpara <- hyper; Dat <- ALL_dat; cur_sam <- ini_sam; N_iter <- N_burn
fn.run.MCMC <- function(hpara, Dat, cur_sam, N_iter)
{
    YY <- Dat$Y
    JJ <- Dat$J
    PP <- Dat$P
    Tot_N <- Dat$Tot_N
    NN <- Dat$n
    Rep_n <- Dat$rep_K
    Miss_ind <- Dat$Miss_ind
    miss_cov_ind <- Dat$miss_cov_ind
    
    MM <- hpara$M
    
    #// prepare for C functions
    YY_c <- array(YY, dim=c(Tot_N*JJ, 1))[,1]
    KK_c <- array(hpara$K, dim=c(Tot_N*hpara$M, 1))[,1]
    
    ### missing will be imputed
    ZZ <- Dat$Z  ## one row for each day
    ZZ_1 <- Dat$Z_1  ## repeated for all replicates - N_Tot rows
    ZZ_1_cov <- ZZ_1[,-c(1:2), drop=FALSE]  ## w/0 index & days - N_Tot rows
    
    print("before burn-in")
    

    for(i_iter in 1:N_iter)
    {
        if((i_iter%%1000)==0)
        {
            print(paste("i.iter=", i_iter))
            print(date())
        }
        
        
        #######################  r_tk ###############################################################
        ### update delta  -- r_tk
        cur_sam$delta <- fn.update.delta(cur_sam$delta, log(cur_sam$psi), cur_sam$w, cur_sam$eta, hpara$c_r, log(cur_sam$r), hpara$u2, Tot_N, hpara$L)
        
        ### update eta  -- r_tk
        cur_sam$eta <- fn.update.eta(cur_sam$delta, hpara$L, cur_sam$w, hpara$u2, hpara$v2_eta, log(cur_sam$r), hpara$c_r)
        
        ### update w  -- r_tk
        cur_sam$w <- fn.update.w(cur_sam$delta, hpara$L, hpara$a_w, hpara$b_w, cur_sam$eta, log(cur_sam$r), hpara$c_r, hpara$u2, cur_sam$w)
        
        ### update psi  -- r_tk
        cur_sam$psi <- fn.update.psi(cur_sam$delta[,1], hpara$d, hpara$L)
        ########################  end of r_tk related parameters ##############################################################################
        
        #######################  th0_j ###############################################################
        ### update delta_th0  -- th0_j
        cur_sam$delta_th0 <- fn.update.delta(cur_sam$delta_th0, log(cur_sam$psi_th0), cur_sam$w_th0, cur_sam$eta_th0, hpara$c_th0, cur_sam$th0, hpara$u2_th0, JJ, hpara$L_th0)
        
        ### update eta_th0  -- th0_j
        cur_sam$eta_th0 <- fn.update.eta(cur_sam$delta_th0, hpara$L_th0, cur_sam$w_th0, hpara$u2_th0, hpara$v2_th0, cur_sam$th0, hpara$c_th0)
        
        ### update w_th0  -- th0_j
        cur_sam$w_th0 <- fn.update.w(cur_sam$delta_th0, hpara$L_th0, hpara$aw_th0, hpara$bw_th0, cur_sam$eta_th0, cur_sam$th0, hpara$c_th0, hpara$u2_th0, cur_sam$w_th0)
        
        ### update psi_th0  -- th0_j
        cur_sam$psi_th0 <- fn.update.psi(cur_sam$delta_th0[,1], hpara$d_th0, hpara$L_th0)
        ######################  end of th0_j related parameters ###############################################################
        
        
        #############################################################################################
        #### update beta, th, th0, s, r
        #############################################################################################
        #### Change Beta  ===> Change BeZ =====> Chage Mu
        #### Change Th  ===> Change KTh =====> Chage Mu
        ### change th0 ===> change mu
        ### update s
        Th_c <- array(cur_sam$th, dim=c(MM*JJ, 1))[,1]
        KTh_c <- array(cur_sam$Kth, dim=c(Tot_N*JJ, 1))[,1]
        
        Beta_c <- array(cur_sam$beta, dim=c(JJ*PP, 1))[,1]
        BeZ_c <- array(cur_sam$beZ, dim=c(Tot_N*JJ, 1))[,1]
        MMu_c <- array(cur_sam$mu, dim=c(Tot_N*JJ, 1))[,1]
        
        Phi_c <- array(cur_sam$phi, dim=c(JJ*PP, 1))[,1]
        
        ZZ_c <- array(ZZ_1_cov, dim=c(Tot_N*PP, 1))[,1]
        
        ### mean of rt=log(r)
        Rt_mean <- c(cur_sam$eta, (hpara$c_r - cur_sam$w*cur_sam$eta)/(1.0-cur_sam$w))
        ind <- cur_sam$delta[,1] + hpara$L*cur_sam$delta[,2]
        Rt_mean <- Rt_mean[ind]
        
        ### mean of th0
        Th0_mean <- c(cur_sam$eta_th0, (hpara$c_th0 - cur_sam$w_th0*cur_sam$eta_th0)/(1.0-cur_sam$w_th0))
        ind <- cur_sam$delta_th0[,1] + hpara$L_th0*cur_sam$delta_th0[,2]
        Th0_mean <- Th0_mean[ind]
        
        output <- .C("fn_update_Beta_Theta", JJ= as.integer(JJ), PP=as.integer(PP), MM=as.integer(MM), Tot_N=as.integer(Tot_N), YY=as.double(YY_c), #
        a_sig=as.double(hpara$a_sig), b_sig=as.double(hpara$b_sig), sig2=as.double(cur_sam$sig2), #
        a_tau=as.double(hpara$a_tau), b_tau=as.double(hpara$b_tau), tau2=as.double(cur_sam$tau2), #
        a_lam=as.double(hpara$a_lam), b_lam=as.double(hpara$b_lam), lam2=as.double(cur_sam$lam2), #
        Rt_m=as.double(Rt_mean), u2=as.double(hpara$u2), RR=as.double(cur_sam$r),
        a_s=as.double(hpara$a_s), b_s = as.double(hpara$b_s), SS=as.double(cur_sam$s), #
        Th0_bar= as.double(Th0_mean), u2_th0=as.double(hpara$u2_th0), Th0=as.double(cur_sam$th0), #
        KK=as.double(KK_c), Th=as.double(Th_c), KTh=as.double(KTh_c),#
        ZZ=as.double(ZZ_c), phi=as.double(Phi_c), Beta=as.double(Beta_c), BeZ=as.double(BeZ_c), #
        MMu=as.double(MMu_c))
        
        ### update th0
        cur_sam$th0 <- output$Th0
        
        ### update th
        cur_sam$th <- array(output$Th, dim=c(MM, JJ))
        cur_sam$Kth <- array(output$KTh, dim=c(Tot_N, JJ))
        
        ### update beta
        cur_sam$beta <- array(output$Beta, dim=c(PP, JJ))   ## p*J matrix
        cur_sam$beZ <- array(output$BeZ, dim=c(Tot_N, JJ))  ## Tot_N*J
        
        cur_sam$mu <- array(output$MMu, dim=c(Tot_N, JJ))  ## Tot_N*J matrix
        
        ### update s
        cur_sam$s <- output$SS
        
        ### update tau2
        cur_sam$tau2 <- output$tau2
        
        ### update sig2
        cur_sam$sig2 <- output$sig2
        
        ### update lam2
        cur_sam$lam2 <- output$lam2
        
        ### update phi
        cur_sam$phi <- array(output$phi, dim=c(PP, JJ))   ## p*J matrix
        
        ### update r
        cur_sam$r <- output$RR
        
        #############################################################################################
        #### END:  update beta, th, th0, s, r
        #############################################################################################
        
        
        ### impute missing Z
        tmp <- fn.impute.Z(YY, NN, Tot_N, JJ, Rep_n, Miss_ind, miss_cov_ind, cur_sam$beta, cur_sam$th0, cur_sam$Kth, cur_sam$r, cur_sam$s, Dat$max_cat, ZZ_1, ZZ, cur_sam$beZ, cur_sam$mu)
        cur_sam$mu <- tmp$MMu
        cur_sam$beZ <- tmp$BeZ
        ZZ <- tmp$Z
        ZZ_1 <- tmp$Z_1
        ZZ_1_cov <- ZZ_1[,-c(1:2),drop=FALSE]   ## w/0 index & days
    }
    
    cur_sam$Z <- ZZ
    cur_sam$Z_1 <- ZZ_1
    cur_sam$Z_1_cov <- ZZ_1_cov

    return(cur_sam)
}


### hpara <- hyper; Dat <- ALL_dat; cur_sam <- ini_sam_1; N_iter <- N_sam
fn.run.MCMC.1 <- function(hpara, Dat, cur_sam, N_iter)
{
    YY <- Dat$Y
    JJ <- Dat$J
    PP <- Dat$P
    Tot_N <- Dat$Tot_N
    NN <- Dat$n
    Rep_n <- Dat$rep_K
    Miss_ind <- Dat$Miss_ind
    miss_cov_ind <- Dat$miss_cov_ind
    miss_sam_ind <- (1:NN)[apply(Miss_ind[,-(1:2),drop=FALSE]==1, 1, sum)>0]  ## miss=1, not miss=0
    N_miss_sam <- length(miss_sam_ind)
    Impute_Z <- array(-1, dim=c(N_miss_sam, length(miss_cov_ind)))
    
    MM <- hpara$M
    
    #// prepare for C functions
    YY_c <- array(YY, dim=c(Tot_N*JJ, 1))[,1]
    KK_c <- array(hpara$K, dim=c(Tot_N*hpara$M, 1))[,1]
    
    
    ### missing will be imputed
    ZZ <- cur_sam$Z ## one row for each day
    ZZ_1 <- cur_sam$Z_1  ## repeated for all replicates - N_Tot rows
    ZZ_1_cov <- cur_sam$Z_1_cov ## w/0 index & days - N_Tot rows
    
    cur_sam$Z <- NULL
    cur_sam$Z_1 <- NULL
    cur_sam$Z_1_cov <- NULL ## w/0 index & days
    
    ####
    SAVE_sam <- NULL
    SAVE_sam$beta <- array(NA, dim=c(PP, JJ, N_iter))
    SAVE_sam$phi <- array(NA, dim=c(PP, JJ, N_iter))
    
    SAVE_sam$s <- array(NA, dim=c(JJ, N_iter))
    SAVE_sam$r <- array(NA, dim=c(Tot_N, N_iter))
    SAVE_sam$th0 <- array(NA, dim=c(JJ, N_iter))
    SAVE_sam$th <- array(NA, dim=c(hpara$M, JJ, N_iter))
    
    SAVE_sam$lam2 <- SAVE_sam$tau2 <- SAVE_sam$sig2 <- array(NA, dim=c(JJ, N_iter))
    
    SAVE_sam$Z_1_cov <- array(NA, dim=c(Tot_N, 10, N_iter))
    
    print("after burn-in")
    
    for(i_iter in 1:N_iter)
    {
        if((i_iter%%1000)==0)
        {
            print(paste("i.iter=", i_iter))
            print(date())
        }

        #######################  r_tk ###############################################################
        ### update delta  -- r_tk
        cur_sam$delta <- fn.update.delta(cur_sam$delta, log(cur_sam$psi), cur_sam$w, cur_sam$eta, hpara$c_r, log(cur_sam$r), hpara$u2, Tot_N, hpara$L)
        
        ### update eta  -- r_tk
        cur_sam$eta <- fn.update.eta(cur_sam$delta, hpara$L, cur_sam$w, hpara$u2, hpara$v2_eta, log(cur_sam$r), hpara$c_r)
        
        ### update w  -- r_tk
        cur_sam$w <- fn.update.w(cur_sam$delta, hpara$L, hpara$a_w, hpara$b_w, cur_sam$eta, log(cur_sam$r), hpara$c_r, hpara$u2, cur_sam$w)
        
        ### update psi  -- r_tk
        cur_sam$psi <- fn.update.psi(cur_sam$delta[,1], hpara$d, hpara$L)
        ########################  end of r_tk related parameters ##############################################################################
        
        #######################  th0_j ###############################################################
        ### update delta_th0  -- th0_j
        cur_sam$delta_th0 <- fn.update.delta(cur_sam$delta_th0, log(cur_sam$psi_th0), cur_sam$w_th0, cur_sam$eta_th0, hpara$c_th0, cur_sam$th0, hpara$u2_th0, JJ, hpara$L_th0)
        
        ### update eta_th0  -- th0_j
        cur_sam$eta_th0 <- fn.update.eta(cur_sam$delta_th0, hpara$L_th0, cur_sam$w_th0, hpara$u2_th0, hpara$v2_th0, cur_sam$th0, hpara$c_th0)
        
        ### update w_th0  -- th0_j
        cur_sam$w_th0 <- fn.update.w(cur_sam$delta_th0, hpara$L_th0, hpara$aw_th0, hpara$bw_th0, cur_sam$eta_th0, cur_sam$th0, hpara$c_th0, hpara$u2_th0, cur_sam$w_th0)
        
        ### update psi_th0  -- th0_j
        cur_sam$psi_th0 <- fn.update.psi(cur_sam$delta_th0[,1], hpara$d_th0, hpara$L_th0)
        ######################  end of th0_j related parameters ###############################################################
        
        
        #############################################################################################
        #### update beta, th, th0, s, r
        #############################################################################################
        #### Change Beta  ===> Change BeZ =====> Chage Mu
        #### Change Th  ===> Change KTh =====> Chage Mu
        ### change th0 ===> change mu
        ### update s
        Th_c <- array(cur_sam$th, dim=c(MM*JJ, 1))[,1]
        KTh_c <- array(cur_sam$Kth, dim=c(Tot_N*JJ, 1))[,1]
        
        Beta_c <- array(cur_sam$beta, dim=c(JJ*PP, 1))[,1]
        BeZ_c <- array(cur_sam$beZ, dim=c(Tot_N*JJ, 1))[,1]
        MMu_c <- array(cur_sam$mu, dim=c(Tot_N*JJ, 1))[,1]
        
        Phi_c <- array(cur_sam$phi, dim=c(JJ*PP, 1))[,1]
        
        ZZ_c <- array(ZZ_1_cov, dim=c(Tot_N*PP, 1))[,1]
        
        ### mean of rt=log(r)
        Rt_mean <- c(cur_sam$eta, (hpara$c_r - cur_sam$w*cur_sam$eta)/(1.0-cur_sam$w))
        ind <- cur_sam$delta[,1] + hpara$L*cur_sam$delta[,2]
        Rt_mean <- Rt_mean[ind]
        
        ### mean of th0
        Th0_mean <- c(cur_sam$eta_th0, (hpara$c_th0 - cur_sam$w_th0*cur_sam$eta_th0)/(1.0-cur_sam$w_th0))
        ind <- cur_sam$delta_th0[,1] + hpara$L_th0*cur_sam$delta_th0[,2]
        Th0_mean <- Th0_mean[ind]
        
        output <- .C("fn_update_Beta_Theta", JJ= as.integer(JJ), PP=as.integer(PP), MM=as.integer(MM), Tot_N=as.integer(Tot_N), YY=as.double(YY_c), #
        a_sig=as.double(hpara$a_sig), b_sig=as.double(hpara$b_sig), sig2=as.double(cur_sam$sig2), #
        a_tau=as.double(hpara$a_tau), b_tau=as.double(hpara$b_tau), tau2=as.double(cur_sam$tau2), #
        a_lam=as.double(hpara$a_lam), b_lam=as.double(hpara$b_lam), lam2=as.double(cur_sam$lam2), #
        Rt_m=as.double(Rt_mean), u2=as.double(hpara$u2), RR=as.double(cur_sam$r),
        a_s=as.double(hpara$a_s), b_s = as.double(hpara$b_s), SS=as.double(cur_sam$s), #
        Th0_bar= as.double(Th0_mean), u2_th0=as.double(hpara$u2_th0), Th0=as.double(cur_sam$th0), #
        KK=as.double(KK_c), Th=as.double(Th_c), KTh=as.double(KTh_c),#
        ZZ=as.double(ZZ_c), phi=as.double(Phi_c), Beta=as.double(Beta_c), BeZ=as.double(BeZ_c), #
        MMu=as.double(MMu_c))
        
        ### update th0
        cur_sam$th0 <- output$Th0
        
        ### update th
        cur_sam$th <- array(output$Th, dim=c(MM, JJ))
        cur_sam$Kth <- array(output$KTh, dim=c(Tot_N, JJ))
        
        ### update beta
        cur_sam$beta <- array(output$Beta, dim=c(PP, JJ))   ## p*J matrix
        cur_sam$beZ <- array(output$BeZ, dim=c(Tot_N, JJ))  ## Tot_N*J
        
        cur_sam$mu <- array(output$MMu, dim=c(Tot_N, JJ))  ## Tot_N*J matrix
        
        ### update s
        cur_sam$s <- output$SS
        
        ### update tau2
        cur_sam$tau2 <- output$tau2
        
        ### update sig2
        cur_sam$sig2 <- output$sig2
        
        ### update lam2
        cur_sam$lam2 <- output$lam2
        
        ### update phi
        cur_sam$phi <- array(output$phi, dim=c(PP, JJ))   ## p*J matrix
        
        ### update r
        cur_sam$r <- output$RR
        
        #############################################################################################
        #### END:  update beta, th, th0, s, r
        #############################################################################################
        
        ### impute missing Z
        tmp <- fn.impute.Z(YY, NN, Tot_N, JJ, Rep_n, Miss_ind, miss_cov_ind, cur_sam$beta, cur_sam$th0, cur_sam$Kth, cur_sam$r, cur_sam$s, Dat$max_cat, ZZ_1, ZZ, cur_sam$beZ, cur_sam$mu)
        cur_sam$mu <- tmp$MMu
        cur_sam$beZ <- tmp$BeZ
        ZZ <- tmp$Z
        ZZ_1 <- tmp$Z_1
        ZZ_1_cov <- ZZ_1[,-c(1:2), drop=FALSE]   ## w/0 index & days
        Impute_Z <- tmp$Impute_Z
        
        
        #### save the current sample
        SAVE_sam$beta[,,i_iter] <- cur_sam$beta
        SAVE_sam$phi[,,i_iter] <- cur_sam$phi
        SAVE_sam$s[,i_iter] <- cur_sam$s
        SAVE_sam$r[,i_iter] <- cur_sam$r
        SAVE_sam$th0[,i_iter] <- cur_sam$th0
        SAVE_sam$th[,,i_iter] <- cur_sam$th
        
        SAVE_sam$lam2[,i_iter] <- cur_sam$lam2
        SAVE_sam$sig2[,i_iter] <- cur_sam$sig2
        SAVE_sam$tau2[,i_iter] <- cur_sam$tau2
        
        SAVE_sam$Z_1_cov[,,i_iter] <- ZZ_1_cov[,1:10]
    }
    
    return(SAVE_sam)
}



### sum_Y: sum of OTU counts per sample
### Tot_N: the number of samples
### std_t: samples
### TT: max time


## YY <- ALL_dat$Y; Tot_N <- ALL_dat$Tot_N; std_t <- ALL_dat$Z_1[,2]; TT <- max(ALL_dat$Z_1[,2]); PP <- ALL_dat$P; JJ <- ALL_dat$J
## YY <- SIM_dat$Y; Tot_N <- SIM_dat$Tot_N; std_t <- SIM_dat$Z_1[,2]; TT <- max(SIM_dat$Z_1[,2]); PP <- SIM_dat$P; JJ <- SIM_dat$J

fn.hyper <- function(YY, Tot_N, std_t, TT, PP, JJ)
{
    ### hyper : fixed hyperparameters
    hpara <- NULL
    
    ### s_j (OTU j): OTU specific dispersion parameter iid Ga(a_s, b_s) with mean a_s/b_s
    hpara$a_s <- 1
    hpara$b_s <- 2
    
    
    #### r_ik: library adjustment.
    #### L : # of mixture components
    hpara$L <- 30
    hpara$u2 <- 1 #0.15  ### each component: w_\ell N(r \eta_\ell, u^2) + (1-w_\ell)N(r, (c_r - w_\ell*\eta_\ell)/(1-w_\ell), u2)
    
    sum_Y <- apply(YY, 1, sum)
    r_mle <- (sum_Y/sum(sum_Y))
    
    hpara$c_r <- mean(log(r_mle))
    
    ### w_\ell \iid \Be(a_w, b_w)
    hpara$a_w <- 1
    hpara$b_w <- 1
    
    ### \eta_\ell \iid \Nor(0, v2_\eta)
    hpara$v2_eta <- 1
    
    #### psi \sim Dir(d, ... , d) -- L-dim
    hpara$d <- 10
    
    #############  end of r_ik
    
    #### th_0j: baseline abundance of otu j.
    #### L_th0 : # of mixture components
    hpara$L_th0 <- 50
    hpara$u2_th0 <- 1  ### each component: w_\ell N(r \eta_\ell, u^2) + (1-w_\ell)N(r, (c_r - w_\ell*\eta_\ell)/(1-w_\ell), u2)
    
    r_mle <- rep.col(r_mle, JJ)
    th0_mle <- log(apply(YY/r_mle, 2, mean))
    
    hpara$c_th0 <- mean(th0_mle)
    
    ### w_\ell \iid \Be(a_w, b_w)
    hpara$aw_th0 <- 1
    hpara$bw_th0 <- 1
    
    ### \eta_\ell \iid \Nor(0, v2_\eta)
    hpara$v2_th0 <- 2.0
    
    #### psi \sim Dir(d, ... , d) -- L-dim
    hpara$d_th0 <- 10
    
    #############  end of th0_j
    
    ############# beta: effects of covariate
    ### beta_jp \iid N(0, \sig2_j \phi_jp) & sig2_j \iid IG(a_sig, b_sig) & phi_jp ~ Exp(\lam2_j/2) & lam2_j ~ Ga(a_lam, b_lam)
    hpara$a_sig <- 3
    hpara$b_sig <- 3
    
    ## mean a/b and variance a/b^2
    hpara$a_lam <- 0.5
    hpara$b_lam <- 0.5
    
    
    ###########  theta: time-dependent part
    hpara$M <- 13 #30  ## # of basis points
    tt0 <- 10
    hpara$u <- seq(-tt0, TT+tt0, length.out=hpara$M)  ## basis points
    hpara$gam2 <-  ((TT+tt0*2)/hpara$M)^2  ## var for the Kernel
    
    ### K: n*M matrix.
    hpara$K <- array(NA, dim=c(Tot_N, hpara$M))
    for(i_m in 1:hpara$M)
    {
        hpara$K[,i_m] <- dnorm(std_t, hpara$u[i_m], sqrt(hpara$gam2))
    }
    
    #plot(hpara$u, hpara$K[1,], type="b")
    #for(i_n in 1:Tot_N)
    #lines(hpara$u, hpara$K[i_n,], type="b")
    
    ### \theta_mj \iid N(0, \tau2_j) & \tau2_j \iid \IG(a_\tau, b_\tau)
    hpara$a_tau <- 0.1
    hpara$b_tau <- 0.1
    
    return(hpara)
}


#### cur_sam
# 1. r: Tot_N-dim vector
# 2. s: J-dim vector
# 3. mu: Tot_N*J matrix:
# 4. delta: Tot_N*2 matrix (indexes for r_tk)
# 5. psi: weights for L components (L-dim vector)
# 6. w: L-dim vector
# 7. eta: L-dim vector
# 8. beta: p*J matrix
# 9. beZ: (p*J)*(Tot_N*p)  ==> Tot_N*J
# 10. th0: J-dim vector
# 11. th: M*J matrix
# 12. Kth: Tot_N*J matrix
# 13: tau0_2: scalar
# 14: tau_2: J-dim vector

# hpara <- hyper; YY <- ALL_dat$Y; ZZ_1 <- ALL_dat$Z_1; Tot_N <- ALL_dat$Tot_N; PP <- ALL_dat$P; JJ <- ALL_dat$J

fn.initialize <- function(hpara, YY, ZZ_1, Tot_N, PP, JJ)
{
    cur_sam <- NULL
    
    ### r
    sum_Y <- apply(YY, 1, sum)
    cur_sam$r <- sum_Y/sum(sum_Y)
    
    ### s_j (OUT j): OTU specific dispersion parameter iid Ga(a_s, b_s)
    cur_sam$s <- rgamma(JJ, hpara$a_s, hpara$b_s)
    
    #### about r_ik
    cur_sam$delta <- array(NA, dim=c(Tot_N, 2))
    cur_sam$delta[,1] <- sample((1:hpara$L), Tot_N, TRUE, prob=rep(1, hpara$L))
    cur_sam$delta[,2] <- sample((0:1), Tot_N, TRUE, prob=c(1, 1))
    
    cur_sam$psi <- rep(1, hpara$L)/hpara$L  ### initally equal prob.
    cur_sam$w <- rbeta(hpara$L, hpara$a_w, hpara$b_w)  #### from the prior
    cur_sam$eta <- rnorm(hpara$L, hpara$c_r, sqrt(hpara$v2_eta))  ### from the prior
    
    ### no effect
    cur_sam$beta <- array(0.01, dim=c(PP, JJ))
    cur_sam$beZ <-  ZZ_1[,-(1:2)]%*%cur_sam$beta # array(0, dim=c(Tot_N, JJ))
    cur_sam$sig2 <- 1/rgamma(JJ, hpara$a_sig, hpara$b_sig) ## var for beta_pj
    cur_sam$lam2 <- rgamma(JJ, hpara$a_lam, hpara$b_lam)
    
    cur_sam$phi <- array(NA, dim=c(PP, JJ))
    for(j in 1:JJ)
    {
        cur_sam$phi[,j] <- rgamma(PP, 1, cur_sam$lam2[j]/2)
    }
    
    
    #### about th_0j
    ### use some empirical values
    r_mat <- rep.col(cur_sam$r, JJ)
    cur_sam$th0 <- log(apply((YY/r_mat), 2, median) + 0.01)
    
    cur_sam$delta_th0 <- array(NA, dim=c(JJ, 2))
    cur_sam$delta_th0[,1] <- sample((1:hpara$L_th0), JJ, TRUE, prob=rep(1, hpara$L_th0))
    cur_sam$delta_th0[,2] <- sample((0:1), JJ, TRUE, prob=c(1, 1))
    
    cur_sam$psi_th0 <- rep(1, hpara$L_th0)/hpara$L_th0  ### initally equal prob.
    cur_sam$w_th0 <- rbeta(hpara$L_th0, hpara$aw_th0, hpara$bw_th0)  #### from the prior
    cur_sam$eta_th0 <- rnorm(hpara$L_th0, hpara$c_th0, sqrt(hpara$v2_th0))  ### from the prior
    
    ### other th
    cur_sam$th <- array(0, dim=c(hpara$M, JJ))
    cur_sam$Kth <- array(0, dim=c(Tot_N, JJ))
    
    th0_mat <- rep.row(cur_sam$th0, Tot_N)
    cur_sam$mu <- exp(th0_mat + cur_sam$Kth + cur_sam$beZ)
    
    cur_sam$tau2 <- 1/rgamma(JJ, hpara$a_tau, hpara$b_tau) ## var for th_mj
    ### \theta_mj \iid N(0, \tau2_j) & \tau2_j \iid \IG(a_\tau, b_\tau)
    
    return(cur_sam)
    
}



################################################################################################################################
##  Delta :  indexes for r_t,k
##  Modeling RRt = log(r_t,k)
##  Tot_N * 2: Column 1: one of L components and Column 2: one for part with eta_\ell and part with (c_r - w_\ell*\eta_\ell)/(1-w_\ell)
################################################################################################################################
################################################################################################################################
##  Delta_th0 :  indexes for th0_j
##  JJ * 2: Column 1: one of L components and Column 2: one for part with eta_\ell and part with (c_r - w_\ell*\eta_\ell)/(1-w_\ell)
################################################################################################################################
# cur_sam$delta <- fn.update.delta(cur_sam$delta, log(cur_sam$psi), cur_sam$w, cur_sam$eta, hpara$c_r, log(cur_sam$r), hpara$u2, Tot_N, hpara$L)
# Delta_cur <- cur_sam$delta; log_Psi <-  log(cur_sam$psi); ww <- cur_sam$w; Eta <- cur_sam$eta; CR <- hpara$c_r; RRt <- log(cur_sam$r); u2 <- hpara$u2; LL <- hpara$L
# Delta_cur <- cur_sam$delta_th0; log_Psi <- log(cur_sam$psi_th0); ww <- cur_sam$w_th0; Eta <- cur_sam$eta_th0; CR <- hpara$c_th0; RRt <- cur_sam$th0; u2 <- hpara$u2_th0; Tot_N <- JJ; LL <- hpara$L_th0
fn.update.delta <- function(Delta_cur, log_Psi, ww, Eta, CR, RRt, u2, Tot_N, LL)
{
    #### Common items
    P1 <- c((log_Psi + log(ww)), (log_Psi + log(1.0-ww)))  ### weights without likelihood
    Mean_tmp <- c(Eta, (CR - ww*Eta)/(1.0-ww))  ### mean for the normal parts
    LL_grid <- (1:(2*LL)) ### first L for (\ell, 0) and the next L for (\ell, 1)
    
    for(i in 1:Tot_N)
    {
        Prob <- P1 - (RRt[i] - Mean_tmp)^2/2.0/u2
        Prob <- exp(Prob - max(Prob))
        
        ind <- sample(LL_grid, 1, FALSE, Prob)
        
        Delta2 <- 1*(ind > LL)
        Delta1 <- ind - (Delta2*LL)
        
        Delta_cur[i, ] <- c(Delta1, Delta2)  ## delta1 \in \{1, ... , L\} and delta2 \in \{0,1\}
    }
    
    return(Delta_cur)
}




################################################################################################################################
## eta_\ell - r_ik
################################################################################################################################
################################################################################################################################
## eta_\ell - th0_j
################################################################################################################################
## Delta <- cur_sam$delta; LL <- hpara$L; u2 <- hpara$u2; v2_eta <- hpara$v2_eta; RRt <- log(cur_sam$r); C_R <- hpara$c_r; ww <- cur_sam$w

fn.update.eta <- function(Delta, LL, ww, u2, v2_eta, RRt, C_R)
{
    w_1_w <- ww/(1-ww)
    Cr_1_w <- C_R/(1.0-ww)
    
    Delta1 <- Delta[,1]  ### delta1 \in \{1, ..., L\}
    Delta2 <- Delta[,2]  ### delta2 \in \{0, 1\}
    
    Eta <- rep(NA, LL)
    
    for(l in 1:LL)
    {
        ind_l0 <- ((Delta1==l)*(Delta2==0))
        ind_l1 <- ((Delta1==l)*(Delta2==1))
        
        v2_tmp <- 1/v2_eta + sum(ind_l0)/u2 + sum(ind_l1)*(w_1_w[l])^2/u2
        v2_tmp <- 1/v2_tmp
        
        m_tmp <- C_R/v2_eta + (sum(ind_l0*RRt) - w_1_w[l]*sum(ind_l1*(RRt - Cr_1_w[l])))/u2
        m_tmp <- v2_tmp*m_tmp
        
        Eta[l] <- rnorm(1, m_tmp, sqrt(v2_tmp))
    }
    
    return(Eta)
}


################################################################################################################################
## w: weights for 0 and 1 for each \ell
################################################################################################################################
###  Delta <- cur_sam$delta; LL <- hpara$L; aa_w <-  hpara$a_w; bb_w <- hpara$b_w; Eta <- cur_sam$eta; C_R <- hpara$c_r; u2 <- hpara$u2; w_cur <- cur_sam$w
fn.update.w <- function(Delta, LL, aa_w, bb_w, Eta, RRt, C_R, u2, w_cur)
{
    Delta1 <- Delta[,1]
    Delta2 <- Delta[,2]
    
    ## propose a new value of w on the logit
    tmp <- log(w_cur/(1.0-w_cur))
    tmp <- tmp + rnorm(LL, 0, 0.5)
    w_pro <- exp(tmp)/(1+exp(tmp))
    
    for(l in 1:LL)
    {
        w_l_cur <- w_cur[l]
        
        ## propose a new value of w on the logit
        w_l_pro <- w_pro[l]
        
        Delta_l0 <- (Delta1==l)*(Delta2==0)
        Delta_l1 <- (Delta1==l)*(Delta2==1)
        RRt_l <- RRt[Delta_l1]
        
        aa_1 <- aa_w + sum(Delta_l0)  ### prior + likelihood + J
        bb_1 <- bb_w + sum(Delta_l1)  ### prior + likelihood + J
        
        m_cur <- (C_R - w_l_cur*Eta[l])/(1-w_l_cur)
        m_pro <- (C_R - w_l_pro*Eta[l])/(1-w_l_pro)
        
        A_cur <- aa_1*log(w_l_cur) + (bb_1)*log(1.0-w_l_cur) - sum((RRt_l - m_cur)^2)/2.0/u2   ## the last term ==0 if nothing is in cluster l
        A_pro <- aa_1*log(w_l_pro) + (bb_1)*log(1.0-w_l_pro) - sum((RRt_l - m_pro)^2)/2.0/u2
        
        if(log(runif(1)) < (A_pro - A_cur))
        {
            w_cur[l] <- w_l_pro
        }## if(log(runif(1)) < (A_pro - A_cur))
    }
    
    
    return(w_cur)
}


################################################################################################################################
## psi: weights for each \ell
################################################################################################################################
# Delta_1 <- cur_sam$delta[,1]; DD <- hpara$d; LL <- hpara$L
fn.update.psi <- function(Delta_1, DD, LL)
{
    ### number of i's in each \ell
    L_cnt <- table(c(Delta_1, (1:LL))) - 1
    
    psi <- rgamma(LL, L_cnt+DD, 1)
    psi <- psi/sum(psi)
    
    return(psi)
}





### repeat x n times in rows ===> produce an n*dim(x) matrix
rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
}

### repeat x n times in columns ===> produce an dim(x)*n matrix
rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


#### K: Tot_N * M where  M: the number of basis points
#### Th0: J-dim vector
#### Th: M*J matrix
####  ==> KTh: Tot_N*J matrix

#### Z: Tot_N*p covariates   *** Note ***: Z is the same for all K_t replicates at time t
## ### \Beta: p*J matrix.
#### BeZ: (p*J)*(Tot_N*p)  ==> Tot_N*J

##> dim(ALL_dat$Z_1)
##[1] 158  12
##[1] 158  14 -- water temp & chlorophyll

#> dim(ALL_dat$Z)
#[1] 56 12


## ZZ_1_cur <- ZZ_1; ZZ_cur <- ZZ; BeZ_cur <- cur_sam$beZ;  MMu_cur <- cur_sam$mu; Beta <- cur_sam$beta; Th0 <- cur_sam$th0; KTh <- cur_sam$Kth; RR <- cur_sam$r; Max_Cat <- Dat$max_cat; SS <- cur_sam$s;
fn.impute.Z <- function(YY, NN, Tot_N, JJ, Rep_n, Miss_ind, miss_cov_ind, Beta, Th0, KTh, RR, SS, Max_Cat, ZZ_1_cur, ZZ_cur, BeZ_cur, MMu_cur)
{
    Th0_Kth <- rep.row(Th0, Tot_N) + KTh ###  Tot_N*J matrix
    RR_mat <- rep.col(RR, JJ)  ## Tot_N*J matrix
    
    ## imputed values
    Z_impute <- array(-1, dim=c(NN, length(miss_cov_ind)))
    
    
    for(i_p in miss_cov_ind)  ##Alex & dino have some missing
    {
        miss_t_i <- (1:NN)[Miss_ind[,(i_p+2)] ==1]  ### missing time points
        
        ### cov indeces
        if(i_p > 1)
        {
            tmp <- sum(Max_Cat[1:(i_p-1)])
        }else{
            tmp <- 0
        }
        cov_ind <- ((tmp+1):(tmp+Max_Cat[i_p]))
        
        ### initialize
        ZZ_cur[miss_t_i, 2+cov_ind] <- 0
        
        ## samples
        for(t_i in miss_t_i)
        {
            ### samples at the time point with missing Z values
            ind_sam <- (1:Tot_N)[ZZ_1_cur[,1]==t_i]       ## indeces from 1 to Tot_N
            
            ### initialize
            ZZ_1_cur[ind_sam, 2+cov_ind] <- 0
            
            Y_i <- YY[ind_sam,]                    # K_i*J matrix
            
            SS_mat <- rep.row(SS, Rep_n[t_i])            ###  K_i*J matrix
            BeZ_mat <- array(NA, dim=c(Rep_n[t_i], JJ, Max_Cat[i_p]+1))
            mu_mat <- array(NA, dim=c(Rep_n[t_i], JJ, Max_Cat[i_p]+1))
            
            ## Z at time point i_t
            ZZ_i <- as.matrix(ZZ_cur[t_i, -(1:2)])  ### cat=0
            Prob_Z <- rep(NA, Max_Cat[i_p]+1)
            
            ### category 0
            i_c <- 0
            ZZ_i_tmp <- ZZ_i
            
            BeZ_mat[,,i_c+1] <- rep.row((ZZ_i_tmp)%*%Beta, Rep_n[t_i])                 # K_i*J matrix
            mu_mat[,,i_c+1] <- mu_tmp <- exp(Th0_Kth[ind_sam,] + BeZ_mat[,,i_c+1])      # K_i*J matrix
            m_tmp <- RR_mat[ind_sam,]*SS_mat*mu_tmp                                     # K_i*J matrix
            
            Prob_Z[i_c+1] <- sum(Y_i*log(m_tmp/(1.0+m_tmp)) - (1/SS_mat)*log(1+m_tmp))
            
            ###  i_c =1, ... max_cat
            for(i_c in 1:Max_Cat[i_p])
            {
                ZZ_i_tmp <- ZZ_i
                ZZ_i_tmp[cov_ind[i_c]] <- 1
                
                BeZ_mat[,,i_c+1] <- rep.row((ZZ_i_tmp)%*%Beta, Rep_n[t_i])                 # K_i*J matrix
                mu_mat[,,i_c+1] <- mu_tmp <- exp(Th0_Kth[ind_sam,] + BeZ_mat[,,i_c+1])      # K_i*J matrix
                m_tmp <- RR_mat[ind_sam,]*SS_mat*mu_tmp                                     # K_i*J matrix
                
                Prob_Z[i_c+1] <- sum(Y_i*log(m_tmp/(1.0+m_tmp)) - (1/SS_mat)*log(1+m_tmp))
            }  ## for(i_c in 1:Max_Cat[i_p])
            
            Prob_Z <- exp(Prob_Z - max(Prob_Z))
            c_sam <- sample((0:Max_Cat[i_p]), 1, FALSE, Prob_Z)
            
            ## if cat > 0
            if(c_sam > 0)
            {
                ZZ_cur[t_i, 2+cov_ind[c_sam]] <- 1
                ZZ_1_cur[ind_sam, 2+cov_ind[c_sam]] <- 1
            } # if(c_sam > 0)
            
            Z_impute[t_i,i_p] <- c_sam
            
            MMu_cur[ind_sam,] <- mu_mat[,,c_sam+1]
            BeZ_cur[ind_sam,] <- BeZ_mat[,,c_sam+1]
            
        } ## for(t_i in miss_t_i)
        
        
    } ## for(i_p in miss_cov_ind)  ##Alex & dino have some missing
    
    return(list(Z_1=ZZ_1_cur, Z=ZZ_cur, BeZ=BeZ_cur, MMu=MMu_cur, Impute_Z = Z_impute))
    
    
}





