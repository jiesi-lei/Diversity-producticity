# calculate rrn copy number for each combination
# 
comm = read.csv("../raw_data/occurrence.csv",row.names = 1)
rrn =  read.csv("data/rrn_tolerance/rrn.csv",row.names = 1)

copy_number <- copy_number.sd <- c()
for (i in 1:81){
rrn_seq = c(t(comm[i,]*rrn$rrn.copy.number))
copy_number = c(copy_number,sum(rrn_seq)/sum(comm[i,]))
copy_number.sd[i] = sd(rrn_seq[rrn_seq > 0])
}

copy_number = data.frame(rownames(comm),copy_number,copy_number.sd)
write.csv(copy_number, "data/rrn_tolerance/copy_number.csv")

#
source("code/0_setup.R")
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
theme_set(theme_test(base_family = "Arial"))

# plot ----
# k ~ species richness
ggplot(res_full) +
  aes(y = k, x =  copy_num, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Maximum'~ OD[600]~'values')) +
  xlab("log(rrn copy number)") +
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
ggsave("plot/k/k_copy_number.pdf", width = 7.25, height = 2.75)

ggplot(filter(res, Species_num == 1)) +
  aes(y = k, x =  log(copy_num), fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Maximum'~ OD[600]~'values')) +
  xlab("log(rrn copy number)") +
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
ggsave("plot/k/k_copy_number_mono.pdf", width = 7.25, height = 2.75)

# r ~ copy number
ggplot(res_full) +
  aes(y = r, x =  copy_num, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Maximum'~ OD[600]~'values')) +
  xlab("log(rrn copy number)") +
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
ggsave("plot/k/k_copy_number.pdf", width = 7.25, height = 2.75)

ggplot(filter(res, Species_num == 1)) +
  aes(y = r, x =  log(copy_num), fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  ylab(bquote('Maximum'~ OD[600]~'values')) +
  xlab("log(rrn copy number)") +
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
ggsave("plot/k/k_copy_number_mono.pdf", width = 7.25, height = 2.75)

fm1 = lm(k ~ log(copy_num), filter(res, trt == "Control"))
summary(fm1)
fm1 = lm(k ~ log(copy_num), filter(res, trt == "Salinity"))
summary(fm1)
fm1 = lm(k ~ I(log(copy_num))^2+log(copy_num), filter(res, trt == "Salinity"))
summary(fm1)
fm1 = lm(k ~ log(copy_num), filter(res, trt == "Starvation"))
summary(fm1)
fm1 = lm(k ~ I(log(copy_num))^2+log(copy_num), filter(res, trt == "Starvation"))
summary(fm1)

ggplot(filter(res, Species_num == "1"), aes(copy_num, r)) + geom_point(aes(color = Species_num)) + 
  geom_smooth(aes(color = Species_num), method = "lm" ) + 
  facet_wrap(vars(trt), ncol = 3, scales = "fixed")

ggplot(filter(res, Species_num == "1"), aes(log(copy_num), k)) + geom_point(aes(color = trt)) + 
  geom_smooth(aes(color = trt), method = "lm") +
  facet_wrap(vars(trt), ncol = 3, scales = "fixed")

ggplot(res, aes(log(copy_num), k)) + geom_point(aes(color = trt)) + 
  geom_smooth(aes(color = trt), method = "lm") +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y")

ggplot(res, aes(copy_num.sd, MPD)) + geom_point(aes(color = trt)) + 
  geom_smooth(aes(color = trt), method = "lm") +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y")

fm1 = lm(MPD ~ log(copy_num), res)
summary(fm1)
