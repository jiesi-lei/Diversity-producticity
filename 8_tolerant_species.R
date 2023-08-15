# Test the effect of tolerant species
rm(list = ls())
source("code/0_data_setup.R")

# Salinity ----
ggplot(salt) +
  aes(y = log_theoretical_r, x =  salt.tol, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Maximal specific growth rate '(min^-1))) +
  xlab("Tol. sp") +
  ylim(0,0.02)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", color = "black", size = 12),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
summary(lm(log_theoretical_r ~ salt.tol, salt))
summary(lm(log_theoretical_r ~ I(salt.tol^2)+salt.tol, salt))

ggplot(res_salt) +
  aes(y = K, x =  salt.tol, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Maximal yield (OD units)')) +
  xlab("Tol. sp") +
  ylim(0,2)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", color = "black", size = 12),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
summary(lm(K ~ salt.tol, salt))


ggplot(salt) +
  aes(y = decrease_rate, x =  salt.tol, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Death rate '(~x10^-4 ~ min^-1))) +
  xlab("Tol. sp") +
  ylim(0,2)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", color = "black", size = 12),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
summary(lm(decrease_rate ~ salt.tol, salt))
summary(lm(decrease_rate ~ I(salt.tol^2)+salt.tol, salt))

# Starvation ----
ggplot(hungry) +
  aes(y = log_theoretical_r, x =  starve.tol, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Maximal specific growth rate '(min^-1))) +
  xlab("Tol. sp") +
  ylim(0,0.04)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", color = "black", size = 12),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
summary(lm(log_theoretical_r ~ starve.tol, hungry))

ggplot(hungry) +
  aes(y = K, x =  starve.tol, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Maximal yield (OD units)')) +
  xlab("Tol. sp") +
  #ylim(0,2)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", color = "black", size = 12),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
summary(lm(K ~ starve.tol, hungry))

ggplot(hungry) +
  aes(y = decrease_rate, x =  starve.tol, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Death rate '(~x10^-4 ~ min^-1))) +
  xlab("Tol. sp") +
  ylim(0,2)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", color = "black", size = 12),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
summary(lm(decrease_rate ~ starve.tol, hungry))

# 只有增加the number of samples才能显著
# species number ----
lm(log_theoretical_r ~ Species_num*salt.tol, res_salt)
fm1 = car::Anova(lm(log_theoretical_r ~ Species_num*salt.tol, res_salt),type = 2)
fm1
lm(K ~ Species_num*salt.tol, res_salt)
fm1 = car::Anova(lm(K ~ Species_num*salt.tol, res_salt),type = 2)
fm1
lm(decrease_rate ~ Species_num*salt.tol, res_salt)
fm1 = car::Anova(lm(decrease_rate ~ Species_num*salt.tol, res_salt),type = 2)
fm1

lm(log_theoretical_r ~ Species_num*salt.tol, res_salt)
fm1 = car::Anova(lm(log_theoretical_r ~ Species_num*salt.tol, res_salt),type = 2)
fm1
lm(K ~ Species_num*salt.tol, res_salt)
fm1 = car::Anova(lm(K ~ Species_num*salt.tol, res_salt),type = 2)
fm1
lm(decrease_rate ~ Species_num*salt.tol, res_salt)
fm1 = car::Anova(lm(decrease_rate ~ Species_num*salt.tol, res_salt),type = 2)
fm1

lm(log_theoretical_r ~ Species_num*starve.tol, res_hungry)
fm1 = car::Anova(lm(log_theoretical_r ~ Species_num*starve.tol, res_hungry),type = 2)
fm1
lm(K ~ Species_num*starve.tol, res_hungry)
fm1 = car::Anova(lm(K ~ Species_num*starve.tol, res_hungry),type = 2)
fm1
lm(decrease_rate ~ Species_num*starve.tol, res_hungry)
fm1 = car::Anova(lm(decrease_rate ~ Species_num*starve.tol, res_hungry),type = 2)
fm1

# species number ----

# linear mixed models
# species number ----
library(lme4)
fm1 = lmer(log_theoretical_r ~ salt.tol + (1|Species_num), res_salt)
fm1
car::Anova(fm1)

fm1 = lmer(K ~ salt.tol + (1|Species_num), res_salt)
fm1
car::Anova(fm1)

fm1 = lmer(decrease_rate ~ salt.tol + (1|Species_num), res_salt)
fm1
car::Anova(fm1)

fm1 = lmer(log_theoretical_r ~ starve.tol + (1|Species_num), res_hungry)
fm1
car::Anova(fm1)

fm1 = lmer(K ~ starve.tol + (1|Species_num), res_hungry)
fm1
car::Anova(fm1)

fm1 = lmer(decrease_rate ~ starve.tol + (1|Species_num), res_hungry)
fm1
car::Anova(fm1)

# PD ----
fm1 = lmer(log_theoretical_r ~ salt.tol + (1|PD), res_salt)
fm1
car::Anova(fm1)

fm1 = lmer(K ~ salt.tol + (1|PD), res_salt)
fm1
car::Anova(fm1)

fm1 = lmer(decrease_rate ~ salt.tol + (1|PD), res_salt)
fm1
car::Anova(fm1)

fm1 = lmer(log_theoretical_r ~ starve.tol + (1|PD), res_hungry)
fm1
car::Anova(fm1)

fm1 = lmer(K ~ starve.tol + (1|PD), res_hungry)
fm1
car::Anova(fm1)

fm1 = lmer(decrease_rate ~ starve.tol + (1|PD), res_hungry)
fm1
car::Anova(fm1)

