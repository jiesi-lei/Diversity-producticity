# coefficients of variation
rm(list = ls())
source("code/0_data_setup.R")
library(dplyr)
library(ggpmisc)

# write.csv(res,"data/res.data.csv")

res = filter(res, decrease_rate < 4 & decrease_rate >= 0)
res$decrease_rate.1 = res$decrease_rate * 10000 #坐标轴上标记*10^4
res.long <- reshape2::melt(res,
                           measure.vars = c("Species_num","PD","MPD"),
                           value.name = "diversity",
                           variable.name = "diversity_type")
# res.long = filter(res.long, decrease_rate > -1e-3)
# 
# res = filter(res, decrease_rate > -1e-3)

my_comparisons <- list(c("Control","Salinity"), c("Salinity", "Starvation"),
                       c("Control", "Starvation"))

ggplot(res.long) +
  aes(y = log_theoretical_r.cv, x =  diversity, fill = trt, colour = trt) + 
  geom_point(shape = "circle", alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2)+x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('CV of maximal specific growth rate')) +
  xlab("Biodiversity") +
  # ylim(0,10)+
  # ggthemes::theme_few() +
  facet_wrap(diversity_type ~ trt, scales = "free", strip.position = "top") +
  # scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  theme_test()+
  theme(strip.text = element_text(size = 5),
        axis.text = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 11, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 11, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1),
        strip.background = element_blank(),
        text = element_text(family = "Arial"))# +
  # stat_poly_eq(formula = y ~ x, 
  #              aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
  #              parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5) +
  # stat_poly_eq(formula = y ~ I(x^2)+x, 
  #              aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
  #              parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5,
  #              vjust = 1.5)
ggsave("plot/cv/growth_rate_cv_diversity.pdf", width = 7.5, height = 7.5)

ggplot(res.long) +
  aes(y = K.cv, x =  diversity, fill = trt, colour = trt) + 
  geom_point(shape = "circle", alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2)+x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('CV of maximal yield')) +
  xlab("Biodiversity") +
  # ylim(0,10)+
  # ggthemes::theme_few() +
  facet_wrap(diversity_type ~ trt, scales = "free", strip.position = "top") +
  # scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  theme_test()+
  theme(strip.text = element_text(size = 5),
        axis.text = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 11, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 11, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1),
        strip.background = element_blank(),
        text = element_text(family = "Arial")) #+
  # stat_poly_eq(formula = y ~ x, 
  #              aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
  #              parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5) 
  # stat_poly_eq(formula = y ~ I(x^2)+x, 
  #              aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
  #              parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5,
  #              vjust = 1.5)
ggsave("plot/cv/yield_cv_diversity.pdf", width = 7.5, height = 7.5)

ggplot(res.long) +
  aes(y = decrease_rate.cv, x =  diversity, fill = trt, colour = trt) + 
  geom_point(shape = "circle", alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2)+x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('CV of death rate')) +
  xlab("Biodiversity") +
  # ylim(0,10)+
  # ggthemes::theme_few() +
  facet_wrap(diversity_type ~ trt, scales = "free", strip.position = "top") +
  # scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  theme_test()+
  theme(strip.text = element_text(size = 5),
        axis.text = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 11, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 11, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1),
        strip.background = element_blank(),
        text = element_text(family = "Arial")) # +
  # stat_poly_eq(formula = y ~ x, 
  #              aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
  #              parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5) +
  # stat_poly_eq(formula = y ~ I(x^2)+x, 
  #              aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
  #              parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5,
  #              vjust = 1.5)
ggsave("plot/cv/death_rate_cv_diversity.pdf", width = 7.5, height = 7.5)

# Faceted - regressions ----
library(lme4) #calculate BIC

# maximal specific growth rate ----
sink("plot/cv/log_theoretical_r.cv.regression.txt")

# species number
fm1 = lm(log_theoretical_r.cv ~ Species_num, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(Species_num^2)+Species_num, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(log_theoretical_r.cv ~ Species_num, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(Species_num^2)+Species_num, filter(res,trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(log_theoretical_r.cv ~ Species_num, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(Species_num^2)+Species_num, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

# phylogenetic diversity
fm1 = lm(log_theoretical_r.cv ~ PD, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(PD^2)+PD, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(log_theoretical_r.cv ~ PD, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(PD^2)+PD, filter(res, trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(log_theoretical_r.cv ~ PD, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(PD^2)+PD, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

# mean pairwise distance
fm1 = lm(log_theoretical_r.cv ~ MPD, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(PD^2)+MPD, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(log_theoretical_r.cv ~ MPD, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(MPD^2)+MPD, filter(res, trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(log_theoretical_r.cv ~ MPD, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(log_theoretical_r.cv ~ I(MPD^2)+MPD, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

closeAllConnections()

# maximal yield ----
sink("plot/cv/maximal_yield.cv.regression.txt")

# species number
fm1 = lm(K.cv ~ Species_num, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(Species_num^2)+Species_num, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K.cv ~ Species_num, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(Species_num^2)+Species_num, filter(res,trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K.cv ~ Species_num, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(Species_num^2)+Species_num, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

# phylogenetic diversity
fm1 = lm(K.cv ~ PD, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(PD^2)+PD, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K.cv ~ PD, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(PD^2)+PD, filter(res, trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K.cv ~ PD, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(PD^2)+PD, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

# mean pairwise distance
fm1 = lm(K.cv ~ MPD, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(PD^2)+MPD, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K.cv ~ MPD, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(MPD^2)+MPD, filter(res, trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K.cv ~ MPD, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K.cv ~ I(MPD^2)+MPD, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

closeAllConnections()

# decrease rate ----
sink("plot/cv/decrease_rate.cv.regression.txt")

# species number
fm1 = lm(decrease_rate.cv ~ Species_num, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(Species_num^2)+Species_num, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(decrease_rate.cv ~ Species_num, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(Species_num^2)+Species_num, filter(res,trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(decrease_rate.cv ~ Species_num, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(Species_num^2)+Species_num, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

# phylogenetic diversity
fm1 = lm(decrease_rate.cv ~ PD, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(PD^2)+PD, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(decrease_rate.cv ~ PD, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(PD^2)+PD, filter(res, trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(decrease_rate.cv ~ PD, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(PD^2)+PD, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

# mean pairwise distance
fm1 = lm(decrease_rate.cv ~ MPD, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(PD^2)+MPD, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(decrease_rate.cv ~ MPD, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(MPD^2)+MPD, filter(res, trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(decrease_rate.cv ~ MPD, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(decrease_rate.cv ~ I(MPD^2)+MPD, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

closeAllConnections()