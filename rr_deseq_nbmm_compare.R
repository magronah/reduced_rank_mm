library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
source("func2.R")
fig_path=   paste0(getwd(), "fig/")
#####################################################
path1   =   paste0(getwd(),"/50_200_previous_500sim/")
path1   =   paste0(getwd(),"/50_200/")
path2   =   paste0(getwd(),"/100_300/")
path3   =   paste0(getwd(),"/150_500/")
#####################################################
dd1     =   load_data(path1) 
dd2     =   load_data(path2) 
dd3     =   load_data(path3) 
#####################################################
pp1     =   do.call(rbind,dd1[["confint"]]) %>% 
              arrange(true_param)
pp1$width =  pp1$upr - pp1$lwr

pp2     =   do.call(rbind,dd2[["confint"]]) %>% 
              arrange(true_param)
pp3     =   do.call(rbind,dd3[["confint"]]) %>% 
              arrange(true_param)

names(pp1)
mod <- c("rr","nbmm")
pp11  = pp1  %>%
        mutate(lwr2 = lwr - average_estimate,
               upr2 = upr - average_estimate)

ggplot(pp11 %>% filter(model %in% mod) , aes(y= param_name, color = model)) +
  #geom_text(aes(x=upr,label=upr,color=model),size=2, hjust=-0.4,vjust=-0.9, show.legend = FALSE)+
  #geom_text(aes(x=lwr,label=lwr,color=model),size=2, hjust=1.2,vjust=-0.9, show.legend = FALSE)+
  geom_errorbarh(aes(xmin=lwr2,xmax=upr2,color=model),height=4)+
  facet_wrap(~model) +
  #geom_point() +
  #geom_line() +
  theme_bw()
###################################################################
# Create plots in a loop
plot_list <- list()

for (i in seq(1, 200, by = 10)) {
  taxon_range <- paste0("taxon", i:(i + 9))
  plot_list[[length(plot_list) + 1]] <- ggplot(pp1 %>% filter(param_name %in% taxon_range), 
                                               aes(x = true_param, y = param_name, color = model)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    coord_flip() +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), position = position_dodge(width = 0.5), height = 0.2) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Estimate", title = paste("Parameters", i, "-", i + 9))
}

(plot_list[[1]] | plot_list[[2]]) / 
(plot_list[[3]] | plot_list[[4]]) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
############################################################
bias1    =    err_extract(dd1, "avg_bias")
bias11   =    reorganise_dd(bias1, "bias")
bias11$label  = factor("50 subjects and 200 ntaxa")

bias2   =   err_extract(dd2, "avg_bias")
bias22   =   reorganise_dd(bias2, "bias")
bias22$label  = factor("100 subjects and 300 ntaxa")

bias3   =   err_extract(dd3, "avg_bias")
bias33   =   reorganise_dd(bias3, "bias")
bias33$label  = factor("150 subjects and 500 ntaxa")
###################################################
mse1    =   err_extract(dd1, "avg_mse")
mse11   =   reorganise_dd(mse1, "mse")
mse11$label  = factor("50 subjects and 200 ntaxa")

mse2   =   err_extract(dd2, "avg_mse")
mse22   =   reorganise_dd(mse2, "mse")
mse22$label  = factor("100 subjects and 300 ntaxa")

mse3   =   err_extract(dd3, "avg_mse")
mse33   =   reorganise_dd(mse3, "mse")
mse33$label  = factor("150 subjects and 500 ntaxa")
###################################################
mse_dd  =  rbind(mse11,mse22,mse33)
bias_dd =  rbind(bias11,bias22,bias33)
##################################################
n = 11
oka_col = c(
  "#0000FF",
  "#556B2F", 
  "#E23D28", 
  "#E69F00", 
  "#000000",
  "#56B4E9", 
  "#D55E00"  # Vermilion
  #"#56B4E9", # Sky Blue
  #"#009E73", # Bluish Green
  #"#F0E442", # Yellow
  #"#0072B2", # Blue
)

mse222 = mse2[mse2$model %in% c("rrzi", "deseq", "zinbmm", "nbmm"), ]
rrzi_value <- mse222$average_value[mse222$model == "rrzi"]

mse222$percentage_change <- ((mse222$average_value - rrzi_value) / rrzi_value) * 100

mse222$model <- factor(mse222$model, levels = c("rrzi", "deseq", "zinbmm", "nbmm"))

g1 <- ggplot(mse222, aes(x = model, y = average_value, color = model)) +
  geom_point(size = 3) +  # Increase point size
  geom_hline(yintercept = rrzi_value, linetype = "dashed", color = "black", 
             linewidth= 1.2) +  # Increase line thickness
  scale_color_manual(values = oka_col) +
  ylab("average value across taxa") +
  custom_theme(n) +
  labs(color = "model") +
  theme(axis.title.x = element_blank())

g1




g1=ggplot(mse33, aes(x =  type, y = average_value, 
                         color = factor(model, levels = unique(model)))) +
  geom_point() +
  scale_color_manual(values = oka_col) +
  ylab("average value across taxa") +
  custom_theme(n) +
  labs(color = "model")+
  #facet_wrap(~label, scales = "free") +
  theme(axis.title.x = element_blank())
g1

dim(dd1$dd$rr)
x11()



g2 =  ggplot(bias_dd, aes(x =  type, y = average_value,  color =model)) +
  geom_point() +
  scale_color_manual(values = oka_col) +
  ylab("average value across taxa") +
  custom_theme(n) +
  labs(color = "model")+
  facet_wrap(~label, scales = "free") +
  theme(axis.title.x = element_blank())
g2


#########################################################################

g3 =  ggplot(mse3, aes(x =  type, y = average_value,  color =model)) +
  geom_point() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())
###############################################################
g4 =  ggplot(bias1, aes(x =  type, y = average_value,  color =model)) +
      geom_point() +
    #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
      ylab("average value across taxa") +
      custom_theme(n) +
      theme(axis.title.x = element_blank())

g5 =  ggplot(bias2, aes(x =  type, y = average_value,  color =model)) +
  geom_point() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())


g6 =  ggplot(bias3, aes(x =  type, y = average_value,  color =model)) +
  geom_point() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())


(g1|g2|g3)/(g4|g5|g6) +   plot_layout(guides = "collect") 
  theme(legend.position = "bottom") # Position the legend


g4

View(dd1$error$rr)

 
g2= ggplot(err1, aes(x =  param_name, y = mse, 
                                    group = model, color = model)) +
  geom_point() +
  geom_line() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange","true_param" = "blue")) +
  xlab("taxa") +
  ylab("mean squared error") +
  #custom_theme(n) +
  theme(axis.text.x = element_blank()) 

g2
##plots


##### comparing average estimate with true parameter value
n  = 13
g1 = ggplot(pp1, aes(x = param_name, y = average_estimate, group = model, 
                    color= factor(model, levels = unique(model) ))) +
  geom_point() +
  geom_point(aes(x = param_name, y = true_param,color="true_param")) +
  geom_line() +
  #scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange", "true_param" = "blue")) +
  labs(color = "model") +   
  xlab("taxa") +
  ylab("average estimate")  +
  #custom_theme(n) +
  theme_bw() +
  theme(axis.text.x = element_blank()) 

ggsave(paste0(fig_path,"deseq_rr_est.png"), plot = g1, width = 8, height = 6, dpi = 500)
###################################################################
names(dd1$error$rr)
##### Compare their MSE
rr_error    =    error_cal(rr, true_param, model = "with_rr")
deseq_error    =    error_cal(deseq_est, true_param, model = "without_rr")

full_summary_dd_all         =    rbind(rr_error$full_summary_dd,    
                                       deseq_error$full_summary_dd)

single_summary_dd_all       =      rbind(rr_error$single_sum_dd,    
                                       deseq_error$single_sum_dd)
###################################################################
g2= ggplot(full_summary_dd_all, aes(x =  param_name, y = mse, 
                                    group = model, color = model)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange","true_param" = "blue")) +
  xlab("taxa") +
  ylab("mean squared error") +
  custom_theme(n) +
  theme(axis.text.x = element_blank()) 

ggsave(paste0(fig_path,"deseq_rr_mse.png"), plot = g2, width = 8, height = 6, dpi = 500)
###################################################################
g3 = ggplot(err1, aes(x =  param_name, y = bias, group = model, color = model)) +
  geom_point() +
  geom_line() +
 # scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange" )) +
  xlab("taxa") +
  ylab("bias") + 
  #custom_theme(n) +
  theme(axis.text.x = element_blank()) 
g3
ggsave(paste0(fig_path,"deseq_rr_bias.png"), plot = g3, width = 8, height = 6, dpi = 500)
###################################################################
g4=ggplot(single_summary_dd_all, aes(x =  type, y = average_value,
                                     group = type,  color =model)) +
  geom_point() +
  scale_color_manual(values = c("with_rr" = "red", "without_rr" = "orange")) +
  ylab("average value across taxa") +
  custom_theme(n) +
  theme(axis.title.x = element_blank())

ggsave(paste0(fig_path,"deseq_rr_avgmse.png"), plot = g4, width = 8, height = 6, dpi = 500)
###################################################################
##confidence intervals
plt  = list(); numbers <- 1:100; group_size <- 10

number_list <- split(numbers, ceiling(seq_along(numbers) / group_size))
position_dodge_width <- 0.5

for(i in 1:length(number_list)){
  elem     =  number_list[[i]]  
  dd_sub   =  dd[dd$param_name %in% paste0("taxon",elem),]
  plt[[i]] =  ggplot(dd_sub, aes(x = param_name, y = true_param, group = model, color = model)) +
    geom_pointrange(aes(ymin = lwr, ymax = upr),
                    position = position_dodge(width = position_dodge_width)) +
    geom_point(position = position_dodge(width = position_dodge_width), color="black",
               size = 2.5) +  
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank()) +  # Remove individual x labels
    xlab("taxa") +
    ylab("true parameter") 
}


# Combine the plots and add a single x-axis label at the bottom
combined_plot <- (plt[[1]] | plt[[2]]) / (plt[[3]] | plt[[4]]) / 
  (plt[[5]] | plt[[6]]) / (plt[[7]] | plt[[8]]) / 
  (plt[[9]] | plt[[10]]) +
  plot_layout(guides = "collect") +
  theme(plot.margin = unit(c(1, 1, 2, 1), "cm"))  # Adjust bottom margin

# Add the x-axis label "taxa" at the bottom
combined_plot + plot_annotation(caption = "taxa") & 
  theme(plot.caption = element_text(hjust = 0.5, vjust = 1, size = 12, face = "bold"))
(plt[[1]]|plt[[2]])/(plt[[3]]|plt[[4]])/(plt[[5]]|plt[[6]])/(plt[[7]]|plt[[8]])/(plt[[9]]|plt[[10]]) +  plot_layout(guides = "collect")

#################################################################################
##coverage 
coverage_cal = function(dd){
  df <- dd %>%
    mutate(
      coverage = if_else(true_param > lwr & true_param < upr, 1, 0)
    )
  df
}

cov = coverage_cal(rr_para_conf)
mean(cov$coverage)
cov2 = coverage_cal(deseq_para_conf)
mean(cov2$coverage)
###################################################################
