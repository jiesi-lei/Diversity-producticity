# Maximal biomass
# Growth rate
source("code/0_data_setup.R")
library(dplyr)
library(ggpmisc)

res.long <- reshape2::melt(res,
                           measure.vars = c("Species_num","PD","MPD"),
                           value.name = "diversity",
                           variable.name = "diversity_type")

# boxplot
my_comparisons <- list(c("Control","Salinity"), c("Salinity", "Starvation"),
                       c("Control", "Starvation"))
ggplot(res) +
  aes(y = K, x =  trt) +
  # geom_boxplot(outlier.shape = NA) +
  geom_violin(aes(fill = trt), trim = FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 1, alpha = 0.5) + 
  geom_boxplot(width = 0.15) +
  # geom_boxplot(width=.1)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = F) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  ylab(bquote('Maximal yield (OD units)')) +
  #ylim(0,0.045)+
  ggthemes::theme_few() +
  xlab("") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
ggsave("plot/maximal_yield/K_violin_tall.pdf", width = 3, height = 5.3)
aggregate(res$K, list(res$trt), FUN=mean)
# 1    Control 1.1797432
# 2   Salinity 1.1108546
# 3 Starvation 0.6443745
(1.1108546/1.1797432-1)*100
(0.6443745/1.1797432-1)*100

# maximal yield ----
ggplot(res.long) +
  aes(y = K, x =  diversity, fill = trt, colour = trt) + 
  geom_point(shape = "circle", alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2)+x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('Maximal yield (OD units)')) +
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
               #aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
               #parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5) +
  #stat_poly_eq(formula = y ~ I(x^2)+x, 
               #aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
               #parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5,
               #vjust = 1.5)
ggsave("plot/maximal_yield/K_diversity.pdf", width = 7.5, height = 7.5)

# Faceted - regressions ----
library(lme4) #calculate BIC
sink("plot/maximal_yield/K_regression.txt")

# species number
fm1 = lm(K ~ Species_num, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(Species_num^2)+Species_num, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K ~ Species_num, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(Species_num^2)+Species_num, filter(res,trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K ~ Species_num, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(Species_num^2)+Species_num, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

# phylogenetic diversity
fm1 = lm(K ~ PD, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(PD^2)+PD, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K ~ PD, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(PD^2)+PD, filter(res, trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K ~ PD, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(PD^2)+PD, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

# mean pairwise distance
fm1 = lm(K ~ MPD, filter(res, trt == "Control"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(MPD^2)+MPD, filter(res, trt == "Control"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K ~ MPD, filter(res, trt == "Salinity"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(MPD^2)+MPD, filter(res, trt == "Salinity"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

fm1 = lm(K ~ MPD, filter(res, trt == "Starvation"))
summary(fm1)
extractAIC(fm1)
fm2 = lm(K ~ I(MPD^2)+MPD, filter(res, trt == "Starvation"))
summary(fm2)
extractAIC(fm2)
extractAIC(fm2)-extractAIC(fm1)
anova(fm1,fm2)

closeAllConnections()
################### End of the new script ####################

