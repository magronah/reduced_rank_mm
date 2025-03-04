library(reformulas)
library(RhpcBLASctl)
library(Matrix)
library(huge)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(glmmTMB)
library(MuMIn)
source("func2.R")
#############################################################
seed = 103; ntaxa = 100; nsubj = 50; beta = 3; betadisp = 0
nsim =  1 ; n_gt = 2;   d    =  2
#############################################################
ntaxa; nsubj; beta; betadisp; 
form <- count ~ 1 + us(1 + group|taxon) + rr(0 + taxon | subject, d=2)
#form <- count ~ 1 + us(1 + group|taxon) + us(0 + taxon | subject)
#####################################################################
dd          =   meta_data(ntaxa, nsubj)
#####################################################################
aicc_list  =   diff =  mod_list  = list()
# tt_eff     =   list(c(1,1,1), c(0.5,0.5,0.5), c(0,0,0.1), c(-1,-1, -0.1))
## glmmTMB interpretes (1,1 as being on the log scale 
## and 1 is 𝜃0/1+𝜃20
## https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html#mappings
#tt_var  =  c(1,1); pho = 0.9; tt_corr   = pho/ sqrt(1 - pho^2)

total_logvar <- 4 #keep change this. This is the only variable thing really
tt_logvar <- c(0.1, 0.1)#c(0.1, 0.1)   
tt_cor <- 0.3       # Correlation term (cor)
rr_logvar <-  c(2,0.5)#c(0.5,2)
set.seed(seed)
#rr_logvar <- rep(2, 4950)#rnorm(4950,mean=-1,sd=1) #c(10, 10)  # rr1 and rr2 (rr_logvar components)
steps <- 10

res <- generate_variance_proportions(total_logvar, tt_logvar,
                                    tt_cor,
                                    rr_logvar, steps)

# Generate results
# res <- generate_increasing_total_variance(tt_logvar,
#                                           tt_cor, rr_logvar,
#                                           steps = steps)

for(i in 1:steps){
  tt_vector <- as.numeric(res[i, c("tt_sd1_cont", "tt_sd2_cont", "tt_cor_cont")])
  rr_vector <-  as.numeric(res[, grep("^rr", names(res))][i,])

  #length(get_theta_corr(ntaxa, nsubj,  seed=seed))
  #theta_true =  c(
   # tt_vector,
    #rr_vector,
    #get_theta_logSD(ntaxa, seed = seed),
    #get_theta_corr(ntaxa, nsubj,  seed=seed)
  #)
  ###################################################################
  theta_true =  c(
    tt_vector,
    get_theta_rr(d, ntaxa,rr_vector,seed =seed)
  )
  true_pars0  =   lst(beta, theta = theta_true, betadisp)
  ####################################################################  
  #### simulate fixed bs
  pars0 <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = seed,
                        newdata    =  dd,
                        newparams  =  true_pars0,
                        family     =  nbinom2, 
                        return_val =  "pars"
  )
  
  stopifnot(theta_true == pars0[names(pars0)=="theta"])
  b_true_all      =   pars0[names(pars0)=="b"]
  true_b_index    =   seq(2,(2*ntaxa),2)
  true_param      =   b_true_all[true_b_index]
  true_b          =   data.frame(param_name = paste0("taxon", 1:ntaxa), 
                                 true_param = true_param)
  ################################################################################
  #now simulate with bs fixed to bval
  true_pars1  =  c(true_pars0, lst(b = b_true_all))
  pars1 <- simulate_new(RHSForm(form, as.form = TRUE), nsim = 1, seed = NULL,
                        newdata    =  dd,
                        newparams  =  true_pars1,
                        family     =  nbinom2, 
                        return_val = "pars"
  )
  stopifnot(b_true_all  == pars1[names(pars1) == "b"])
  ####################################################################  
  count_list <- simulate_new(RHSForm(form, as.form = TRUE),
                             newdata   =  dd,
                             newparams =  true_pars1,
                             family    =  nbinom2,
                             nsim      =  1,
                             seed   =  seed)[[1]]
  #####################################################################
  dd$count  =   count_list
  df        =   otu_meta_fun(dd)
  #####################################################################
  gprior  <- data.frame(prior = "gamma(2, 2.5)",
                        class = "theta_sd",
                        coef = "")
  #####################################################################
  form2   <- count ~ 1 + us(1 + group|taxon) + rr(0 + taxon | subject, d=2) + offset(normalizer)
  #form2  =   update(form,.~. + offset(normalizer))
  form3  =   count ~ 1 + us(1 + group|taxon) + offset(normalizer)
  #####################################################################
  options(glmmTMB_openmp_debug = TRUE)
  blas_set_num_threads(1)
  
  par_ctrl <- glmmTMBControl(
    parallel = list(n = 10, autopar = TRUE)
  )
  #####################################################################
  mod1  <- tryCatch({
    glmmTMB(form2, data = df,
            family  = nbinom2,
            prior   = gprior,
            REML    = FALSE,
            control = par_ctrl)
  }, error = function(e) {
    message("Error in first attempt, trying again without prior...")
    glmmTMB(form2, data = df,
            family  = nbinom2,
            REML    = FALSE,
            control = par_ctrl)
  })
  #####################################################################
  mod2 = glmmTMB(form3, data = df,
                 family  = nbinom2,
                 prior   = gprior,
                 REML    = FALSE,
                 control = par_ctrl)
  #####################################################################
  aicc_list[[i]] = c(rr = AICc(mod1), us =   AICc(mod2))
  diff[[i]]           =  AICc(mod1) -   AICc(mod2)
  mod_list[[i]]       =  list(rr = mod1, us = mod2)
                         ##mod1 = rr; mod2 =us
}

######################################################################
ddf <- data.frame(aicc =unlist(diff))
size = 3; nn= 14; width =  10; height =  5; dpi = 300
ddf$sim  =   factor(paste0("sim",1:steps),
                    levels = paste0("sim",1:steps))

prr = ggplot(ddf, aes(x = sim, y= aicc, group = 1)) +
  geom_point(size = size) +
  geom_line() +
  labs(x =" ", y = "AICc(RR) - AICc(US)") +
  custom_theme(nn) 

ggsave("fig/aic_diff.png", plot = prr, 
       width = width, 
       height = height, 
       dpi = dpi)

diff
res
ds
######################################################################
## power, for each have to do at least 300 simulations and find the proportion
## of them that has pval < 0.05

ds=data.frame(ttvar = apply(res, 1, function(x){
  x["tt_sd1_cont"] + x["tt_sd2_cont"] + x["tt_cor_cont"] }), 
  rrvar = apply(res, 1, function(x){
    x["rr1_sd_cont"] +  x["rr2_sd_cont"]})
) 


ds2 = data.frame(type = rep(c("ttvar","rrvar"), each=steps),  
                 proportion = c(apply(res, 1, function(x){
                   x["prop_tt"]}), 
                   apply(res, 1, function(x){
                     x["prop_rr"]})),
                 sim  =   factor(rep(paste0("sim",1:steps),2),
                                 levels = paste0("sim",1:steps))
) 

ds1  =  data.frame(ttvar = apply(res, 1, function(x){
  x["prop_tt"]}), 
  rrvar = apply(res, 1, function(x){
    x["prop_rr"]})) 


size = 3; nn= 14; width =  10; height =  5; dpi = 300
pr = ggplot(ds2, aes(x = sim, y= proportion, group = type, color =  type)) +
  geom_point(size = size) +
  geom_line() +
  labs(x =" ") +
  scale_color_manual(name = "Total variance proportion", 
                     values = c("ttvar" = "red", "rrvar" = "blue"), 
                     labels = c("ttvar" = "treatment variance",
                                "rrvar" = "reduced rank variance")) +
  custom_theme(nn) 

ggsave("fig/var_prop.png", plot = pr, 
       width = width, 
       height = height, 
       dpi = dpi)
