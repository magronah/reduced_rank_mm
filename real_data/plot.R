library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(AICcmodavg)
##############################################################
path1   =   paste0(getwd(),"/real_data/")
source(paste0(path1,"/fun.R"))
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
######################################################################
# Load models for autism data
autism_path <- paste0(getwd(), "/real_data/autism_data/results/")
atlass_path <- paste0(getwd(), "/real_data/atlass_data/results/")
crohn_path <- paste0(getwd(), "/real_data/CrohnD_data/results/")
soil_path   <- paste0(getwd(), "/real_data/soil_data/results/")

# Load models for autism data
autism_models <- load_models(autism_path, filenames)
atlass_models <- load_models(atlass_path, filenames)
crohn_models <- load_models(crohn_path, filenames)
soil_models   <- load_models(soil_path, filenames)

# Assigning names 
mod_names <- c("RR","RRzi","US","USzi","Nbmm","Zinbmm","DE")
names(autism_models)    =   names(atlass_models)  =    mod_names
names(crohn_models)    =   names(soil_models)    =    mod_names
mod_list    =   lst(autism_models,atlass_models,crohn_models,soil_models)
##################################################################################
res     =  lapply(mod_list, function(x){
                           pp  =  list(RR = x$RR, RRzi = x$RRzi,
                                       US = x$US, USzi =x$USzi)
                           
                           dd   <-  data.frame(aictab(pp, modnames = names(pp)))
                           df1  <-  dd %>% select(Modnames,AICc)
                           na_indices <- which(is.na(df1$AICc))
                           
                           if (length(na_indices) > 0) {
                             df1$AICc[na_indices] <- sapply(na_indices, function(i) my_aicc_fun(pp[[df1$Modnames[i]]]))
                           }
                           df2  <- data.frame(Modnames =  c("Nbmm", "Zinbmm", "DE"),
                                              AICc     =  c(Nbmm   =  x$Nbmm,
                                                            Zinbmm =  x$Zinbmm,
                                                              DE   =   x$DE))
                           df   =   rbind(df1, df2)
                           df = bind_rows(df1, df2) %>%
                             mutate(Delta_AICc = AICc - min(AICc)) %>%
                             arrange(Delta_AICc)
                           df
                           })

#' the 4th data set fitting zi for each taxa leads to worse AICc
#' Nothing changes for the 3rd dataset in terms of the order but each zi is better
#' 
#############################################################
plt   =   list()
nam   =    c("Autism data",
             "Atlass data", 
             "Crohn disease data", 
             "Soil data")

for(i in 1:length(res)){
  ddf      =      res[[i]]
  ddf$Delta_AICc[ddf$Delta_AICc == 0] <- 1e-6
  plt[[i]] = ggplot(ddf, aes(Modnames,Delta_AICc)) + 
            geom_point() +
             scale_y_log10() +
              custom_theme(12) +
             labs(x ="", y = "Change in AICc", title = nam[i])
}

(plt[[1]]|plt[[2]])/(plt[[3]]|plt[[4]])

#############################################################
# Load models for autism data
filename        <-  "nbmm_mod.rds"
autism_nbmm_mod <-  load_models(autism_path, filename)
atlass_nbmm_mod <-  load_models(atlass_path, filename)
crohn_nbmm_mod <-  load_models(crohn_path, filename)
soil_nbmm_mod   <-  load_models(soil_path, filename)

ntaxa_vec    <-    c(autism_mod  =  length(autism_nbmm_mod[[1]]$responses),
                   atlass_mod  =   length(atlass_nbmm_mod[[1]]$responses),
                   crohn_mod  =  length(crohn_nbmm_mod[[1]]$responses),
                   soil_mod    =  length(soil_nbmm_mod[[1]]$responses))

mean_count_list <- list(
                  autism_mod  =  readRDS(paste0(autism_path,"mean_count.rds")),
                  atlass_mod  =   readRDS(paste0(atlass_path,"mean_count.rds")),
                  crohn_mod  =   readRDS(paste0(crohn_path,"mean_count.rds")),
                  soil_mod    =   readRDS(paste0(soil_path,"mean_count.rds"))
                  )

View(dd_full)
dd_full   =   wald_confint(mod = crohn_models$RR)

grp_dd    <-  extract_grp_effect(dd_full,ntaxa)

ggplot(grp_dd,aes(param_name, width, group=1)) +
  geom_point() +
  geom_line() +
  custom_theme(12)


# Add model type
dd$model =   rep("wald",nrow(dd))


library(zCompositions)
library(compositions)
library(corrplot)
library(reshape2)
library(ggplot2)
library(dplyr)
#############################################################
# Step 1: Create Synthetic Microbiome Dataset
set.seed(123)
# Parameters
n_sites <- 5
n_groups <- 2
n_taxa <- 10
n_samples_per_site <- 20

# Create site, group, and sample IDs
site <- rep(paste0("Site", 1:n_sites), each = n_samples_per_site)
group <- rep(rep(c("Control", "Treatment"), each = n_samples_per_site / 2), times = n_sites)
sample_id <- paste0("Sample", 1:(n_sites * n_samples_per_site))

# Generate taxon counts with site-specific patterns
taxa_counts <- matrix(0, nrow = length(site), ncol = n_taxa)
colnames(taxa_counts) <- paste0("Taxon", 1:n_taxa)

for (taxon in 1:n_taxa) {
  for (s in 1:n_sites) {
    site_samples <- which(site == paste0("Site", s))
    mean_count <- 5 + s * 2 + rnorm(1, 0, 2)
    taxa_counts[site_samples, taxon] <- rpois(length(site_samples), lambda = mean_count)
  }
}

# Introduce group effects (e.g., treatment increases certain taxa counts)
group_effect_taxa <- c(1, 2)  # Group effects on Taxon1 and Taxon2
for (taxon in group_effect_taxa) {
  taxa_counts[group == "Treatment", taxon] <- taxa_counts[group == "Treatment", taxon] +
    rpois(sum(group == "Treatment"), lambda = 5)
}

# Create a data frame
microbiome_data <- as.data.frame(taxa_counts)
microbiome_data$site <- site
microbiome_data$group <- group
microbiome_data$sample_id <- sample_id

# Step 2: Apply CLR Transformation to Address Compositionality
# Replace zeros using Bayesian multiplicative replacement
taxa_counts_no_zeros <- cmultRepl(taxa_counts, label = 0, method = "CZM")

# CLR transformation
clr_transformed <- apply(taxa_counts_no_zeros, 1, function(row) clr(as.numeric(row)))
clr_transformed <- t(clr_transformed)  # Transpose to match original structure

# Convert CLR-transformed data to a data frame
clr_data <- as.data.frame(clr_transformed)
colnames(clr_data) <- paste0("Taxon", 1:n_taxa)
clr_data$site <- microbiome_data$site
clr_data$group <- microbiome_data$group
clr_data$sample_id <- microbiome_data$sample_id

# Step 3: Explore Taxon Correlations Across Sites (CLR-Transformed Data)
correlation_by_site_clr <- lapply(unique(clr_data$site), function(s) {
  site_data <- clr_data[clr_data$site == s, 1:n_taxa]
  cor(site_data)
})
View(clr_data)
#############################################################
crohn_path <- paste0(getwd(), "/real_data/CrohnD_data/results/")

crohn_dir <- dirname(crohn_path)
source(paste0(crohn_dir,"/prep_data.R"))
#############################################################
# Load required packages
library(zCompositions)
library(compositions)
library(corrplot)
library(reshape2)
library(ggplot2)
library(dplyr)
library(glmmTMB)

# Ensure your data contains columns: count, age, group, taxon, subject, normalizer
microbiome_data        <-  as.matrix(t(pp$data$countdata))
microbiome_data        <-   cbind(microbiome_data, pp$data$meta_data)
stopifnot(rownames(microbiome_data) == pp$data$meta_data$subject)
rownames(microbiome_data)  <- NULL

# Step 2: CLR Transformation to Address Compositionality
# Replace zeros using Bayesian multiplicative replacement
taxa_counts_no_zeros <- cmultRepl((microbiome_data[, grep("taxon", colnames(microbiome_data))]), 
                                  label = 0, method = "CZM")

# CLR transformation
clr_transformed <- apply(taxa_counts_no_zeros, 1, function(row) clr(as.numeric(row)))
clr_transformed <- t(clr_transformed)  # Transpose to match original structure

# Convert CLR-transformed data to a data frame
clr_data <- as.data.frame(clr_transformed)
colnames(clr_data) <- grep("taxon", colnames(microbiome_data), value = TRUE)
clr_data$age <- microbiome_data$age
clr_data$group <- microbiome_data$group
clr_data$subject <- microbiome_data$subject
clr_data$normalizer <- microbiome_data$normalizer

subject_data = cor(clr_data[,grep("taxon", colnames(clr_data))])

# Step 5: Visualize Correlations Across Subjects (CLR-Transformed Data)
correlation_by_subject <- lapply(unique(clr_data$subject), function(s) {
  subject_data <- clr_data[clr_data$subject == s, grep("taxon", colnames(clr_data))]
  cor(subject_data)
})



# Visualize correlation matrices
par(mfrow = c(2, 3))
for (i in seq_along(correlation_by_subject)) {
  corrplot(correlation_by_subject[[i]], method = "color", 
           title = paste0("CLR-Subject ", i), mar = c(0, 0, 2, 0))
}

# Step 6: Check for Group-Level Effects (CLR-Transformed Data)
# Melt CLR-transformed data for group-level comparison
clr_melted_data <- melt(clr_data, id.vars = c("age", "group", "subject"), 
                        variable.name = "taxon", value.name = "clr_count")

# Calculate mean CLR-transformed abundance by group
group_means_clr <- clr_melted_data %>%
  group_by(group, taxon) %>%
  summarise(mean_clr = mean(clr_count), .groups = "drop")

# Visualize group differences for each taxon
ggplot(group_means_clr, aes(x = taxon, y = mean_clr, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean CLR-Transformed Taxon Abundance by Group", 
       x = "Taxon", y = "Mean CLR Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


















# Step 2: Apply CLR Transformation to Address Compositionality
# Replace zeros using Bayesian multiplicative replacement
dd1     =     dd
rownames(dd1) <- NULL

taxa_count1 =  as.matrix(dd1)
taxa_counts_no_zeros <- cmultRepl(taxa_count1, label = 0, method = "CZM")
#############################################################
# CLR transformation
clr_transformed <- apply(taxa_counts_no_zeros, 1, function(row) clr(as.numeric(row)))
clr_transformed <- t(clr_transformed)  # Transpose to match original structure

# Convert CLR-transformed data to a data frame
clr_data <- as.data.frame(clr_transformed)
colnames(clr_data) <- paste0("Taxon", 1:n_taxa)
clr_data$site <- microbiome_data$site
clr_data$group <- microbiome_data$group
clr_data$sample_id <- microbiome_data$sample_id
#############################################################

# Step 3: Explore Taxon Correlations Across Sites (CLR-Transformed Data)
correlation_by_site_clr <- lapply(unique(clr_data$site), function(s) {
  site_data <- clr_data[clr_data$site == s, 1:n_taxa]
  cor(site_data)
})

# Visualize correlation matrices
par(mfrow = c(2, 3))
for (i in seq_along(correlation_by_site_clr)) {
  corrplot(correlation_by_site_clr[[i]], method = "color", 
           title = paste0("CLR-Site ", i), mar = c(0, 0, 2, 0))
}

# Step 4: Check for Group-Level Effects (CLR-Transformed Data)
# Melt CLR-transformed data for group-level comparison
clr_melted_data <- melt(clr_data, id.vars = c("site", "group", "sample_id"), 
                        variable.name = "taxon", value.name = "clr_count")

# Calculate mean CLR-transformed abundance by group
group_means_clr <- clr_melted_data %>%
  group_by(group, taxon) %>%
  summarise(mean_clr = mean(clr_count), .groups = "drop")

# Visualize group differences for each taxon
ggplot(group_means_clr, aes(x = taxon, y = mean_clr, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean CLR-Transformed Taxon Abundance by Group", 
       x = "Taxon", y = "Mean CLR Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
