library(tibble)
library(dplyr)
library(scam)
library(tidyverse)
library(here)
library(gratia)
library(patchwork)

source("reproducible/func.R")
############################################
nsubj = 150; ntaxa = 500 
path <- paste0(getwd(), "/", nsubj, "_", ntaxa, "/")
##################################################
## Read all the needed data
true_param =  readRDS(paste0(path,"true_param.rds"))

# Read mean count data
mean_count <- readRDS(paste0(path, "mean_count.rds"))

## Read estimates from fitting ensembles
dd   = list(
  rr     =  readRDS(paste0(path, "rr.rds")),
  rrzi   =  readRDS(paste0(path, "rrzi.rds")),
  us     =  readRDS(paste0(path, "us.rds")),
  uszi   =  readRDS(paste0(path, "uszi.rds"))
)
##########################################################
effect_size <-  true_param$true_param

# Compute p-values for all models
pvals <- lapply(dd, pvalue_cal)
names(pvals) <- c("rr", "rrzi", "us", "uszi")
##########################################################  
# Fit GAM models for each p-value set
mod_rr    <- gam_fit(pvals[["rr"]], effect_size, 
                     mean_count, alpha_level = 0.05)  # no error

shrinkfun <- function(x, eps = 0.001) (x + eps)/(1+2*eps)
dd <- tibble(mean_count, effect_size, pval = shrinkfun(pvals[["rr"]]),
             pval_rrzi = shrinkfun(pvals[["rrzi"]])) |>
    mutate(pval_reject = as.numeric(pval < 0.05),
           pval_rrzi_reject = as.numeric(pval_rrzi < 0.05))


brkvec <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.99)

gg0 <- ggplot(dd, aes(x=mean_count, y = abs(effect_size))) +
    scale_x_log10()

gg1 <- gg0 + geom_point(aes(colour = pval)) +
        scale_colour_viridis_c(trans = "logit", breaks = brkvec)
gg2 <- gg0 + geom_point(aes(colour = factor(pval_reject))) +
    scale_colour_manual(values = c("red", "black"))

## since pval_rrzi_reject doesn't give us a problem with scam, let's look at
## which values are *different* ...
with(dd, table(pval_reject, pval_rrzi_reject))

gg3 <- gg0 + geom_point(aes(colour = factor(pval_rrzi_reject),
                            size = factor(pval_reject != pval_rrzi_reject))) +
    scale_colour_manual(values = c("red", "black")) +
    scale_size_manual(values = c(1, 5))

## fit GAM - results look reasonable?
myfit = mgcv::gam(pval_reject ~ te(log(mean_count), abs(effect_size), bs="cr"),
                  data = dd, family = binomial)
gratia:draw(myfit)

## tedmi = double monotone increasing -- seems plausible

if (FALSE) {
    
    myfit = scam::scam(pval_reject ~ s(log(mean_count), abs(effect_size), bs="tedmi"),
                       data = dd, family = binomial, control = list(trace = TRUE, print.warn = TRUE,
                                                                    devtol.fit = 1e-4,
                                                                    steptol.fit = 1e-4,
                                                                    maxit = 1000))
    ## fails at computing eigen(S.t); S.t is completely full of NaN values at this point
    ## regularize?

    ## try various perturbed versions of the data
    
    set.seed(101)
    dd_boot <- dd[sample(nrow(dd)),]
    myfit = scam::scam(pval_reject ~ s(log(mean_count), abs(effect_size), bs="tedmi"),
                       data = dd_boot, family = binomial, control = list(trace = TRUE, print.warn = TRUE,
                                                                         devtol.fit = 1e-4,
                                                                         steptol.fit = 1e-4,
                                                                         maxit = 1000))

    set.seed(101)
    dd_jitter <- dd |> mutate(across(c(mean_count, effect_size, pval), ~ jitter(.x, factor = 1)),
                              pval_reject = as.numeric(pval < 0.05))
    myfit = scam::scam(pval_reject ~ s(log(mean_count), abs(effect_size), bs="tedmi"),
                       data = dd_jitter, family = binomial, control = list(trace = TRUE, print.warn = TRUE,
                                                                           devtol.fit = 1e-4,
                                                                           steptol.fit = 1e-4,
                                                                           maxit = 1000))
}

mod_rrzi  <- gam_fit(pvals[["rrzi"]], effect_size, 
                     mean_count, alpha_level = 0.05)  # no error

mod_us    <- gam_fit(pvals[["us"]], effect_size, 
                     mean_count, alpha_level = 0.05)   # no error

mod_uszi  <- gam_fit(pvals[["uszi"]], effect_size, 
                     mean_count, alpha_level = 0.05)


#######################################################
any(is.na(pvals[["uszi"]]))
any(is.infinite(pvals[["uszi"]]))
#######################################################
#For nsubj = 150; ntaxa = 500, I get the same error for all the model 
## except for mod_us which works fine. 
