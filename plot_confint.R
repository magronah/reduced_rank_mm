library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(here)
source("func2.R")
fig_path=   "fig/"
########################################################
path    =   "200_600/confint/"
titles  =   "100 subjects per group and 600 taxa"
########################################################
# Define filenames to load
filenames <- c(
  "rr.rds",
  "rrzi.rds",
  "us.rds",
  "uszi.rds",
  "nbmm.rds",
  "zinbmm.rds",
  "deseq.rds"  
)

# Load models for autism data
dd        <-  load_models(path, filenames)
names(dd)  =    sub("\\.rds$", "", filenames)
names(dd)  <- c("RR","RRzi","US","USzi","NB","ZNB","DE" )
#####################################################################
nrw   =   600;  ncl = 6
dd_nbmm     <-  discard(dd[["NB"]],   ~ any(dim(.x) != c(nrw, ncl)))
dd_zinbmm   <-  discard(dd[["ZNB"]], ~ any(dim(.x) != c(nrw, ncl)))
######################################################################
common_names1  <-   intersect(names(dd_nbmm), names(dd_zinbmm))
common_names   <-   intersect(common_names1, names(dd$RRzi))

dd_nbmm       <-   dd_nbmm[common_names]
dd_zinbmm     <-   dd_zinbmm[common_names]
dd_deseq      <-   dd$DE[common_names]
dd_rr         <-   dd$RR[common_names]
sub_dd        <-   dd[grep("US|RR", names(dd))]
us_rr         <-   lapply(sub_dd, function(x){x[common_names]})
dd_list       <-   list(NB  =  dd_nbmm, ZNB = dd_zinbmm)
######################################################################
## compute the average proportion of significant taxa
sig1   <-   lapply(us_rr, function(y){x <- map_dfc(y, ~ .x$pvalue)
                                    apply(x, 2, p.adjust, method = "BH")})

sig2   <-   lapply(dd_list, function(y){x <- map_dfc(y, ~ .x$pvalue)
                           apply(x, 2, p.adjust, method = "BH")})

sig3   <-   list(DE = map_dfc(dd_deseq, ~ pull(.x, padj)))

sig   <-    c(sig1, sig2, sig3)
######################################################################
alpha  = 0.05
ddf = data.frame(lapply(sig, function(x){mean(colMeans(x < alpha))}))
ddf =  t(ddf) %>%
       data.frame() %>%
       setNames("prop") %>%
       rownames_to_column("model")

##############################################################
size = 3; nn= 13; width =  7; height =  5; dpi = 300
plt0 = ggplot(ddf, aes(model, prop)) +
  geom_point() +
  custom_theme(nn) +
  labs(x =" ", y = "Average number of p-values less than 0.05")

ggsave("fig/pvalhits.png", plot = plt0, 
       width = width, 
       height = height, 
       dpi = dpi)
##############################################################
##confidence width and coverage plot
true_param    =   readRDS("200_600/true_param.rds")
##################################################################
cov_list  <-   c(dd_list, us_rr, DE = list(dd_deseq))
##################################################################
ddd <- lapply(cov_list, function(x){
  ddf <- lapply(x, function(y){
    df         =   left_join(true_param,y, by = "param_name")
    df$cov     =   ifelse(df$lwr  < df$true_param & df$true_param < df$upr,1, 0)
    df
  })
  names(ddf)  = names(x)
  ddf
})
##############################################################
covv   <-   lapply(ddd, function(y){map_dfc(y, ~ .x$cov)})
dd_cov   <-   data.frame(lapply(covv, function(x){mean(rowMeans(x))}))
ddf1 =  t(dd_cov) %>%
   as.data.frame() %>%
  setNames("coverage") %>%
  rownames_to_column("model")

plt1 = ggplot(ddf1, aes(model, coverage)) +
  geom_point() +
  custom_theme(nn) +
  labs(x =" ", y = "coverage") +
  ggtitle("Average coverage for all taxa")


ggsave("fig/coverage.png", plot = plt1, 
       width = width, 
       height = height, 
       dpi = dpi)

#########################################################################
width_dfs <- lapply(cov_list, function(model) {
  width_df <- as.data.frame(lapply(model, function(sim) pull(sim, width)))
  wdmean   <-  mean(rowMeans(width_df))
  colnames(width_df) <- names(model)
  return(wdmean)
})
width_dd   <- t(data.frame(width_dfs))  %>%
   as.data.frame() %>%
  setNames("width") %>%
  rownames_to_column("model")

plt2 = ggplot(width_dd, aes(model, width)) +
  geom_point() +
  custom_theme(nn) +
  ggtitle("Average confidence width across simulations for all taxa") +
  labs(x =" ", y = "Confidence width")


ggsave("fig/cov_width.png", plot = plt2, 
       width = width, 
       height = height, 
       dpi = dpi)
##################################################################
pl=plt2|plt1
ggsave("fig/cov_nwidth.png", plot = pl, 
       width = 13, 
       height = height, 
       dpi = dpi)

