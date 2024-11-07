#setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
path        =   paste0(getwd(),"/reproducible/new_sim/real_data/")
##############################################################
## have nodes =  1 and 
# cpu-per-task= more
library(tidyr)
library(tidyverse)
library(tibble)
library(DESeq2)
library(RhpcBLASctl)
source("load_glmmTMB.R")
#library(glmmTMB)/
source(paste0(path,"fun.R"))
source(paste0(path,"autism/load_data.R"))
##############################################################
dd_long     =   df_long(dd_filt)
dd_long     =   left_join(dd_long, met_dd, by ="subject")  
###########################################################
dds        =   DESeqDataSetFromMatrix(dd_filt,met_dd, ~group)
dds        =   DESeq(dds,sfType ="poscounts",minReplicatesForReplace=Inf) 
normalizer =   data.frame(normalizer = sizeFactors(dds)) %>% 
               rownames_to_column("subject")

df      =    left_join(dd_long, normalizer, by ="subject")  
###########################################################
form      =   count ~ 1 + us(1 + group|taxon) +  
                     rr(0 + taxon | subject,2) + 
                     offset(normalizer)

##Fit the model
par_ctrl <- glmmTMBControl(
  parallel = list(n = 6, autopar = TRUE),
  optArgs  = list(maxit = 10000)
)



gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")

options(glmmTMB_openmp_debug = TRUE)
blas_set_num_threads(1)

tt = system.time(
    fit1  <-  glmmTMB(form, data = df,
                   family  = nbinom2, 
                   prior   = gprior,
                   REML    = TRUE,
                   control = par_ctrl
  )
)
saveRDS(tt, file=paste0(path,"autism/data/rr_runtime",".rds"))
saveRDS(fit, file=paste0(path,"autism/data/rr_mod",".rds"))



tt = system.time(
  fit <- glmmTMB(formula = form, 
                 data = df, 
                 family = nbinom2, 
                 ziformula = ~1 + (1 | taxon),
                 prior = gprior, 
                 REML = TRUE, 
                 control = par_ctrl)
)


betazi    =  fit$fit$par[["betazi"]]
betadisp  =  
theta     =  fit$fit$par[names(fit$fit$par) == "theta"]

unique(names(fit$fit$par))
#################################
# Generate data
ntaxa    =  max(dim(dd_filt))
grp_pars =  list(name_vec =  as.vector(unique(met_dd$group)), 
                size_vec  =  c(10,10))

param    =   list(ntaxa = 100, 
                  beta_pars =  list(beta= fixef(fit),
                                    betadisp = fit$fit$par[["betadisp"]]),
                  theta     =  theta[names(fit$fit$par) == "theta"],
                  grp_pars  =  list(name_vec = c("control", "treat"), 
                                  size_vec = c(10,10)))

                  

meta_dd <- met_data(ntaxa, grp_pars)

# Simulate new data using the specified parameters
sim_count = simulate_new(
  RHSForm(form, as.form = TRUE),
  newdata =   meta_dd,
  newparams = pars,
  family = nbinom2,
  ziformula = ~ taxon,
  nsim = nsim,
  seed = seed
)


names(fit$fit$par)
#next try the other datasets as well 
#coef(fit)
# library(glmmTMB)
# setting OpenMP: threads = 6L, autopar = TRUE
# setting OpenMP: threads = 1L, autopar = FALSE
# user   system  elapsed 
# 2070.020    6.624  391.032

# setting OpenMP: threads = 6L, autopar = TRUE
# setting OpenMP: threads = 1L, autopar = FALSE
# user   system  elapsed 
# 1743.617    6.990  331.135 




