seed = 101; ntaxa = 600; nsubj = 200; beta = 3; betadisp = 0 
#betazi = -2
nsim =  500 ; n_gt = 2

theta_true =  c(
  get_theta_logSD(n_gt, seed = seed),
  get_theta_corr(n_gt, n_gt,  seed=seed),
  get_theta_logSD(ntaxa, seed = seed),
  get_theta_corr(ntaxa, nsubj,  seed=seed)
)

