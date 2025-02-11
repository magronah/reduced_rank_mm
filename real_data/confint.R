library(tidyverse)
library(foreach)
library(doParallel)
#################################################
path1   =   paste0(getwd(),"/real_data/")
source(paste0(path1,"/fun.R"))
#################################################
# Load models for autism data
autism_path <- paste0(getwd(), "/real_data/autism_data/results/")
atlass_path <- paste0(getwd(), "/real_data/atlass_data/results/")
crohn_path <- paste0(getwd(), "/real_data/CrohnD_data/results/")
soil_path   <- paste0(getwd(), "/real_data/soil_data/results/")
#########################################################
# Define filenames to load
filenames <- c(
  "rr_mod.rds",
  "rrzi_each_mod.rds",
  "us_mod.rds",
  "uszi_each_mod.rds",
  "nbmm_aicc.rds",
  "zinbmm_aicc.rds",
  "deseq_aicc.rds"  
)
#################################################
# Load models for autism data
autism_models <- load_models(autism_path, filenames)
atlass_models <- load_models(atlass_path, filenames)
crohn_models <- load_models(crohn_path, filenames)
soil_models   <- load_models(soil_path, filenames)

# Assigning names 
mod_names <- c("RR","RRzi","US","USzi","Nbmm","Zinbmm","DE")
names(autism_models)    =   names(atlass_models)  =    mod_names
names(crohn_models)    =   names(soil_models)    =    mod_names
mod_list    =   lst(autism_models,atlass_models,
                    crohn_models, soil_models)
###########################################################
us_rr_models <- lapply(mod_list, function(x){
  x[grep("^(RR|US)", names(x))]
})
###########################################################
autism_us_rr  =   us_rr_models$autism_models
atlass_us_rr  =   us_rr_models$atlass_models
crohn_us_rr   =   us_rr_models$crohn_models
 soil_us_rr   =   us_rr_models$soil_models
############################################################
registerDoParallel(cores = max(min(detectCores()-2, 10 ),1)) 
test  =  list(autism_us_rr)

time_taken <- system.time({
res <- lapply(test, function(models) {
   foreach(mod = models, .combine = list) %dopar% wald_confint(mod)
 })
})

print(time_taken)
#is there a way to speed it up?
###########################################################

res <- foreach::foreach(i = 1:length(autism_us_rr)) %dopar% {
  mod   =   autism_us_rr[[i]]
  wald_confint(mod)
}

res <- foreach::foreach(i = 1:length(atlass_us_rr)) %dopar% {
   mod   =   atlass_us_rr[[i]]
   wald_confint(mod)
}

names(us_rr_models)
#################################################

confint_dd_list <- list()
for(i in 1:length(us_rr_models)){
  mod_list    <-   us_rr_models[[i]]
  confint_dd_list[[i]]  <- lapply(mods, function(x){
    wald_confint(mod = x)})
}


