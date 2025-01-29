library(glmmTMB)
library(ggplot2)
library(dplyr)
library(AICcmodavg)
##############################################################
path   =   paste0(getwd(),"/real_data/soil_data/results/")
##############################################################
     rr       =   readRDS(paste0(path,"rr_mod.rds"))
     rrzi     =   readRDS(paste0(path,"rrzi_mod.rds"))
     us       =   readRDS(paste0(path,"us_mod.rds"))
     uszi     =   readRDS(paste0(path,"uszi_mod.rds"))
     nbmm     =   readRDS(paste0(path,"nbmm_aicc.rds"))
   zinbmm     =   readRDS(paste0(path,"zinbmm_aicc.rds"))
    deseq     =   readRDS(paste0(path,"zinbmm_aicc.rds"))
##############################################################
    mod     =   lst(rr,rrzi,us,uszi)
      dd    =   data.frame(aictab(mod, modnames = names(mod)))
##############################################################
df  <- dd %>%
        select(Modnames,AICc, Delta_AICc)
      
a   =  min(df$AICc) -  zinbmm
if(z < 0){
  df 
}
  
###############################################################
load_models <- function(path, filenames) {
  file_paths <- paste0(path, filenames)
  lapply(file_paths, readRDS)
}

# Define filenames to load
filenames <- c(
  "rr_mod.rds",
  "rrzi_mod.rds",
  "us_mod.rds",
  "uszi_mod.rds",
  "nbmm_aicc.rds",
  "zinbmm_aicc.rds",
  "deseq_aicc.rds"  
)

# Load models for soil data
soil_path <- paste0(getwd(), "/real_data/soil_data/results/")
soil_models <- load_models(soil_path, filenames)

# Load models for autism data
autism_path <- paste0(getwd(), "/real_data/autism_data/results/")
autism_models <- load_models(autism_path, filenames)

# Assigning meaningful names for the loaded models
names(soil_models) <- c("rr", "rrzi", "us", "uszi", "nbmm", "zinbmm", "deseq")
names(autism_models) <- c("rr", "rrzi", "us", "uszi", "nbmm", "zinbmm", "deseq")
##################################################################################




# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)


cd = c("autism dataset", "soil dataset", "chrons disease dataset", "atlas dataset")

aic_values <- data.frame(
  Dataset = paste("Dataset", 1:5),
  Model1 = runif(5, 100, 150),
  Model2 = runif(5, 100, 150),
  Model3 = runif(5, 100, 150),
  Model4 = runif(5, 100, 150),
  Model5 = runif(5, 100, 150)
)

# Calculate differences in AIC relative to the best model for each dataset
aic_diff <- aic_values %>%
  pivot_longer(cols = starts_with("Model"), names_to = "Model", values_to = "AIC") %>%
  group_by(Dataset) %>%
  mutate(AIC_Diff = AIC - min(AIC))

# Visualize differences in AIC
ggplot(aic_diff, aes(x = Dataset, y = AIC_Diff, color = Model, group = Model)) +
  #geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Differences in AIC for Models Across Datasets",
    x = "Dataset",
    y = "Difference in AIC",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
      
      
            