library(lattice)
library(nlme)
library(BayesRegMicrobiome)
load("SIM_dat.RData")
#######################################################
# SIM_dat$Y
hyper <- fn.hyper(SIM_dat$Y, SIM_dat$Tot_N, SIM_dat$Z_1[,2], 
                  
                  max(SIM_dat$Z_1[,2]), SIM_dat$P, SIM_dat$J)

names(hyper)
