# 9 . Stability
source("code/0_data_setup.R")
library(dplyr)
library(ggpubr)
library(ggpmisc)

# Get unique Species_num values
species_nums <- c(1,2,4,8,12,16)

# For growth rate -----
# Initialize an empty dataframe to store the differences
result_dataframe_r <- data.frame()

# Iterate over each Species_num value
for (species_num in species_nums[1:6]) {
  # Get the log_theoretical_r values for control, salinity, and hungry samples at the current Species_num
  control_r <- res_normal$log_theoretical_r[res_normal$Species_num == species_num]
  salinity_r <- res_salt$log_theoretical_r[res_salt$Species_num == species_num]
  hungry_r <- res_hungry$log_theoretical_r[res_hungry$Species_num == species_num]
  
  # Perform t-test
  t_result_1 <- t.test(control_r, salinity_r)
  mean_control <- t_result_1$estimate[[1]]
  mean_salinity <- t_result_1$estimate[[2]]
  diff_lb_1 <- t_result_1$conf.int[[1]]
  diff_ub_1 <- t_result_1$conf.int[[2]]
  t_value_1 <- t_result_1$statistic[[1]]
  p_value_1 <- t_result_1$p.value[[1]]
  
  t_result_2 <- t.test(control_r, hungry_r)
  mean_hungry <- t_result_2$estimate[[2]]
  diff_lb_2 <- t_result_2$conf.int[[1]]
  diff_ub_2 <- t_result_2$conf.int[[2]]
  t_value_2 <- t_result_2$statistic[[1]]
  p_value_2 <- t_result_2$p.value[[1]]
  
  # Append the results to the result_dataframe
  result_dataframe_r <- rbind(result_dataframe_r, data.frame(Species_num = c(species_num), 
                                                         mean_control = mean_control,
                                                         mean_salinity = mean_salinity,
                                                         diff_salinity_control = 100*(mean_salinity-mean_control)/mean_control,
                                                         diff_lb_1 = diff_lb_1,
                                                         diff_ub_1 = diff_ub_1,
                                                         diff_salinity_control_ub = -100*diff_lb_1/mean_control,
                                                         diff_salinity_control_lb = -100*diff_ub_1/mean_control,
                                                         t_value_1 = t_value_1, 
                                                         p_value_1 = p_value_1,
                                                         mean_hungry = mean_hungry,
                                                         diff_hungry_control = 100*(mean_hungry-mean_control)/mean_control,
                                                         diff_lb_2 = diff_lb_2,
                                                         diff_ub_2 = diff_ub_2,
                                                         diff_hungry_control_ub = -100*diff_lb_2/mean_control,
                                                         diff_hungry_control_lb = -100*diff_ub_2/mean_control,
                                                         t_value_2 = t_value_2, 
                                                         p_value_2 = p_value_2))
}

# Print the resulting dataframe
print(result_dataframe_r)
write.csv(result_dataframe_r,"results_R/stability/result_dataframe_r.csv")

result_r <- reshape2::melt(result_dataframe_r,
                           measure.vars = c("diff_salinity_control","diff_hungry_control"),
                           value.name = "diversity",
                           variable.name = "diversity_type")

# For yield -----
# Initialize an empty dataframe to store the differences
result_dataframe_yield <- data.frame()

# Iterate over each Species_num value
for (species_num in species_nums[1:6]) {
  # Get the K values for control, salinity, and hungry samples at the current Species_num
  control_r <- res_normal$K[res_normal$Species_num == species_num]
  salinity_r <- res_salt$K[res_salt$Species_num == species_num]
  hungry_r <- res_hungry$K[res_hungry$Species_num == species_num]
  
  # Perform t-test
  t_result_1 <- t.test(control_r, salinity_r)
  mean_control <- t_result_1$estimate[[1]]
  mean_salinity <- t_result_1$estimate[[2]]
  diff_lb_1 <- t_result_1$conf.int[[1]]
  diff_ub_1 <- t_result_1$conf.int[[2]]
  t_value_1 <- t_result_1$statistic[[1]]
  p_value_1 <- t_result_1$p.value[[1]]
  
  t_result_2 <- t.test(control_r, hungry_r)
  mean_hungry <- t_result_2$estimate[[2]]
  diff_lb_2 <- t_result_2$conf.int[[1]]
  diff_ub_2 <- t_result_2$conf.int[[2]]
  t_value_2 <- t_result_2$statistic[[1]]
  p_value_2 <- t_result_2$p.value[[1]]
  
  # Append the results to the result_dataframe
  result_dataframe_yield <- rbind(result_dataframe_yield, data.frame(Species_num = c(species_num), 
                                                                     mean_control = mean_control,
                                                                     mean_salinity = mean_salinity,
                                                                     diff_salinity_control = 100*(mean_salinity-mean_control)/mean_control,
                                                                     diff_lb_1 = diff_lb_1,
                                                                     diff_ub_1 = diff_ub_1,
                                                                     diff_salinity_control_ub = -100*diff_lb_1/mean_control,
                                                                     diff_salinity_control_lb = -100*diff_ub_1/mean_control,
                                                                     t_value_1 = t_value_1, 
                                                                     p_value_1 = p_value_1,
                                                                     mean_hungry = mean_hungry,
                                                                     diff_hungry_control = 100*(mean_hungry-mean_control)/mean_control,
                                                                     diff_lb_2 = diff_lb_2,
                                                                     diff_ub_2 = diff_ub_2,
                                                                     diff_hungry_control_ub = -100*diff_lb_2/mean_control,
                                                                     diff_hungry_control_lb = -100*diff_ub_2/mean_control,
                                                                     t_value_2 = t_value_2, 
                                                                     p_value_2 = p_value_2))
}

# Print the resulting dataframe
print(result_dataframe_yield)

print(result_dataframe_yield)
write.csv(result_dataframe_yield,"results_R/stability/result_dataframe_yield.csv")

# Plot ----
res_growth <- read.csv("results_R/stability/growth_rate_long.csv")
res_yield <- read.csv("results_R/stability/yield_long.csv")

ccodes = c("#6aa84f","#61a0d9")

library(ggplot2)
library(ggpubr)
library(lme4)
library(egg)
library(RColorBrewer)
library(MASS)
library(sfsmisc)
library(dplyr)
library(ggsci)
library("scales")
library(wesanderson)
library(ggh4x)#set scales for facet_wrap
theme_set(ggthemes::theme_few(base_family = "Arial"))

res_growth$Treatment = factor(res_growth$Treatment, levels = c("Salinity","Hungry"))
res_growth$p_value_sig <- (res_growth$p_value < 0.05)
res_growth$Species_num <- as.numeric(res_growth$Species_num)

ggplot(res_growth, aes(x=Species_num, y=relative_change, col=Treatment, shape = p_value_sig)) +#, shape=sig, alpha=sig
  geom_abline(intercept=0, slope=0, col="black", linetype = "dashed") +
  geom_line(aes(group=Treatment)) +
  scale_shape_manual(values=c(1,16)) +
  geom_errorbar(aes(ymin=change_lb, ymax=change_ub),
                width=0) +
  ylab(" ") +
  xlab("Species richness") +
  ylim(-75,50) +
  geom_point(size=3) + 
  scale_colour_manual(values=ccodes) +
  ggthemes::theme_few() +
  theme(axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(color = "black", size = 12, family = "Arial"),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, family = "Arial", margin = margin(t = 8, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 0.75),
        legend.key.size = unit(0.75, 'lines'),
        legend.text = element_text(size=12, family = "Arial"),
        legend.spacing.y = unit(0.75, 'lines'))
ggsave("plot/stability/growth_rate_with_sig_barplot.pdf", width = 4.3, height = 2.9)

res_yield$Treatment = factor(res_yield$Treatment, levels = c("Salinity","Hungry"))
res_yield$p_value_sig <- (res_yield$p_value < 0.05)
res_yield$Species_num <- as.numeric(res_yield$Species_num)

ggplot(res_yield, aes(x=Species_num, y=relative_change, col=Treatment, shape = p_value_sig)) +#, shape=sig, alpha=sig
  geom_abline(intercept=0, slope=0, col="black", linetype = "dashed") +
  geom_line(aes(group=Treatment)) +
  scale_shape_manual(values=c(1,16)) +
  geom_errorbar(aes(ymin=change_lb, ymax=change_ub),
                width=0) +
  ylab(" ") +
  xlab("Species richness") +
  ylim(-75,50) +
  geom_point(size=3) + 
  scale_colour_manual(values=ccodes) +
  ggthemes::theme_few() +
  theme(axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(color = "black", size = 12, family = "Arial"),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, family = "Arial", margin = margin(t = 8, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 0.75),
        legend.key.size = unit(0.75, 'lines'),
        legend.text = element_text(size=12, family = "Arial"),
        legend.spacing.y = unit(0.75, 'lines'))
ggsave("plot/stability/yield_with_sig_barplot.pdf", width = 4.3, height = 2.9)

