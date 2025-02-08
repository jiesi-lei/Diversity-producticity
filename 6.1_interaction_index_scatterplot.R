# interaction index
rm(list = ls())
source("code/0_data_setup.R")

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
library(viridis)
theme_set(ggthemes::theme_few(base_family = "Arial"))

comm = read.csv("../raw_data/occurrence.csv",row.names = 1)

# growth rate ----
filter(res, Species_num == 1&trt == "Control") %>% dplyr::select(.,log_theoretical_r) -> mono.ctrl
id = which(res_normal$Species_num_comb %in% rownames(comm))
comm.ctrl = merge(res_normal$Species_num_comb, comm, )
exptd.r.ctrl = apply(comm,1, function(x){
  sum = sum(x*mono.ctrl)
  count = sum(x)
  exp.r = sum/count
  exp.r
})

filter(res, Species_num == 1&trt == "Salinity") %>% dplyr::select(.,log_theoretical_r) -> mono.salt 
exptd.r.salt = apply(comm,1, function(x){
  sum = sum(x*mono.salt)
  count = sum(x)
  exp.r = sum/count
  exp.r
})

filter(res, Species_num == 1&trt == "Starvation") %>% dplyr::select(.,log_theoretical_r) -> mono.stv
exptd.r.stv = apply(comm,1, function(x){
  sum = sum(x*mono.stv)
  count = sum(x)
  exp.r = sum/count
  exp.r
})

res$exptd.r = c(exptd.r.ctrl,exptd.r.salt,exptd.r.stv)
res$interaction.r = res$log_theoretical_r/res$exptd.r

# subset species richness > 1
res %>% filter(Species_num >1) %>% group_by(trt, Species_num) %>% 
  summarise(interaction.index.r = mean(interaction.r), ii.sd = sd(interaction.r))
res %>% filter(Species_num >1) %>% group_by(trt) %>% 
  summarise(interaction.index.r = mean(interaction.r), ii.se = std.error(interaction.r)) -> b


# plot
ggplot(filter(res, Species_num >1)) + 
  geom_point(aes(x = exptd.r, y = log_theoretical_r, color = Species_num, shape = trt), size =1, stroke = 0.4) +
  scale_shape_manual(values = c(0,1,2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
  scale_color_viridis(option = "D") +
  xlab(bquote("Expected growth rate "(min^-1)))+
  ylab(bquote("Community growth rate "(min^-1)))+
  ggthemes::theme_few()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
res %>% filter(Species_num >1) %>% group_by(trt, Species_num) %>% 
  summarise(interaction.index.r = mean(interaction.r), ii.sd = sd(interaction.r))ggsave("plot/interaction/expected.vs.empirical.r.pdf", width = 4.4, height = 2.9)

mean(res$log_theoretical_r)

# maximal yield ----
filter(res, Species_num == 1&trt == "Control") %>% dplyr::select(.,K) -> mono.ctrl 
comm = read.csv("../raw_data/occurrence.csv",row.names = 1)
id = which(res_normal$Species_num_comb %in% rownames(comm))
comm.ctrl = merge(res_normal$Species_num_comb, comm, )
exptd.k.ctrl = apply(comm,1, function(x){
  sum = sum(x*mono.ctrl)
  count = sum(x)
  exp.k = sum/count
  exp.k
})

filter(res, Species_num == 1&trt == "Salinity") %>% dplyr::select(.,K) -> mono.salt 
exptd.k.salt = apply(comm,1, function(x){
  sum = sum(x*mono.salt)
  count = sum(x)
  exp.k = sum/count
  exp.k
})

filter(res, Species_num == 1&trt == "Starvation") %>% dplyr::select(.,K) -> mono.stv
exptd.k.stv = apply(comm,1, function(x){
  sum = sum(x*mono.stv)
  count = sum(x)
  exp.k = sum/count
  exp.k
})

res$exptd.k = c(exptd.k.ctrl,exptd.k.salt,exptd.k.stv)
res$interaction.k = res$K/res$exptd.k

# subset species richness > 1
res %>% filter(Species_num >1) %>% group_by(trt, Species_num) %>% summarise(interaction.index.k = mean(interaction.k), ii.sd = sd(interaction.k))
res %>% filter(Species_num >1) %>% group_by(trt) %>% 
  summarise(interaction.index.k = mean(interaction.k), ii.se = std.error(interaction.k)) -> b

ggplot(filter(res, Species_num >1)) + 
  geom_point(aes(x = exptd.k, y = K, color = Species_num, shape = trt), size = 1, stroke = 0.4) +
  scale_shape_manual(values = c(0,1,2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
  scale_color_viridis(option = "D") +
  xlab("Expected yield (OD units)")+
  ylab("Community yield (OD units)")+
  ggthemes::theme_few()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
ggsave("plot/interaction/expected.vs.empirical.k.pdf", width = 4.21, height = 2.7)

# death rate ----
filter(res, Species_num == 1&trt == "Control") %>% dplyr::select(.,decrease_rate) -> mono.ctrl 
comm = read.csv("../raw_data/occurrence.csv",row.names = 1)
id = which(res_normal$Species_num_comb %in% rownames(comm))
comm.ctrl = merge(res_normal$Species_num_comb, comm, )
exptd.d.ctrl = apply(comm,1, function(x){
  sum = sum(x*mono.ctrl)
  count = sum(x)
  exp.d = sum/count
  exp.d
})

filter(res, Species_num == 1&trt == "Salinity") %>% dplyr::select(.,decrease_rate) -> mono.salt 
exptd.d.salt = apply(comm,1, function(x){
  sum = sum(x*mono.salt)
  count = sum(x)
  exp.d = sum/count
  exp.d
})

filter(res, Species_num == 1&trt == "Starvation") %>% dplyr::select(.,decrease_rate) -> mono.stv
exptd.d.stv = apply(comm,1, function(x){
  sum = sum(x*mono.stv)
  count = sum(x)
  exp.d = sum/count
  exp.d
})

res$exptd.d = c(exptd.d.ctrl,exptd.d.salt,exptd.d.stv)
res$interaction.d = res$decrease_rate/res$exptd.d

# subset species richness > 1
res %>% filter(Species_num >1) %>% group_by(trt, Species_num) %>% summarise(interaction.index.d = mean(interaction.d), ii.sd = sd(interaction.d))

ggplot(filter(res, Species_num >1)) + 
  geom_point(aes(x = exptd.d, y = decrease_rate, color = Species_num, shape = trt), size = 1, stroke = 0.4) +
  scale_shape_manual(values = c(0,1,2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
  scale_color_viridis(option = "D") +
  xlab(bquote("Expected decrease rate "(min^-1)))+
  ylab(bquote("Community decrease rate "(min^-1)))+
  ggthemes::theme_few()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
ggsave("plot/interaction/expected.vs.empirical.d.pdf", width = 4.21, height = 2.7)

#### End of the new script ####

ggplot(filter(res, Species_num == 2)) + 
  geom_point(aes(x = exptd.k, y = K, color = Combination, shape = trt)) +
  scale_shape_manual(values = c(0,1,2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
  scale_color_viridis(option = "C", discrete = T) +
  xlab("Expected growth rate")+
  ylab("Community growth rate")+
  ggthemes::theme_few()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

filter(res, Species_num == 2)%>% group_by(trt, Combination) %>% 
  summarise(interaction.index.r = mean(interaction.r), ii.sd = sd(interaction.r))-> a

filter(res, Species_num > 1)%>% group_by(trt, Species_num_comb) %>% 
  summarise(interaction.index.k = mean(interaction.k), ii.sd = sd(interaction.k), 
            interaction.index.r = mean(interaction.r), ii.sd = sd(interaction.r), MPD = mean(MPD),
            rrn = mean(copy_num), rrn.sd = mean(copy_num.sd),
            richness = mean(Species_num )) %>% filter(interaction.index.r < 3) -> b

cor.test(b$richness, b$rrn.sd)
cor.test(b$richness, b$rrn)

c = filter(b, trt == "Control")
d = filter(b, trt == "Salinity")
e = filter(b, trt == "Starvation")

cor.test(c$interaction.index.k, c$MPD, method = "pearson")
cor.test(c$interaction.index.r, c$MPD, method = "pearson")
cor.test(c$interaction.index.r, c$rrn, method = "pearson")
cor.test(c$interaction.index.k, c$rrn, method = "pearson")
cor.test(c$interaction.index.k, c$rrn.sd)
cor.test(c$interaction.index.r, c$rrn.sd)
cor.test(d$interaction.index.k, d$MPD, method = "pearson")
cor.test(d$interaction.index.r, d$MPD, method = "pearson")
cor.test(d$interaction.index.r, d$rrn)
cor.test(d$interaction.index.k, d$rrn)
cor.test(d$interaction.index.k, d$rrn.sd)
cor.test(d$interaction.index.r, d$rrn.sd)
cor.test(e$interaction.index.k, e$MPD, method = "pearson")
cor.test(e$interaction.index.r, e$MPD, method = "pearson")
cor.test(e$interaction.index.r, e$rrn)
cor.test(e$interaction.index.k, e$rrn)
cor.test(e$interaction.index.k, e$rrn.sd)
cor.test(e$interaction.index.r, e$rrn.sd)
cor.test(e$interaction.index.r, e$richness)
cor.test(e$interaction.index.k, e$richness)
cor.test(e$MPD, e$rrn.sd)
cor.test(c$MPD, c$rrn.sd)

ggplot(b, aes(x = rrn.sd, y = interaction.index.r)) + 
  geom_point(shape = "circle", alpha = 1, size = 2, color = "gray50") +
  xlab("sd of rrn copy number") +
  ylab("Interaction index for growth rate")+
  geom_smooth(method = "lm", formula = y ~ x, size = 0.75, se = F, color = "black") +
  #scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylim(0.5,2.5) +
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "fixed") +
  # theme_test()+
  theme(strip.text = element_text(size = 12, family = "Arial"),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(color = "black", size = 12, family = "Arial", margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, family = "Arial", margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
ggsave("plot/interaction/r_rrn_sd.pdf", width = 7.25, height = 3)

ggplot(b, aes(x = rrn.sd, y = interaction.index.k)) + 
  geom_point(shape = "circle", alpha = 1, size = 2, color = "gray50") +
  xlab("sd of rrn copy number")+
  ylab("Interaction index for community yield")+
  geom_smooth(method = "lm", formula = y ~ x, size = 0.75, se = F, color = "black") +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylim(0.5,2)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "fixed") +
  # theme_test()+
  theme(strip.text = element_text(size = 12, family = "Arial"),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(color = "black", size = 12, family = "Arial", margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, family = "Arial", margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))
ggsave("plot/interaction/k_rrn_sd.pdf", width = 7.25, height = 3)

summary(lm(interaction.index.k ~ trt*rrn.sd*MPD, b))

summary(lm(interaction.index.r ~ rrn.sd*MPD*trt, b))

library(ggpubr)
my_comparisons <- list(c("Control","Salinity"), c("Salinity", "Starvation"),
                       c("Control", "Starvation"))
ggplot(b) +
  aes(y = interaction.index.k, x =  trt, colour = trt) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2.5, alpha = 0.7) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  ylab(bquote('Specific growth rate '(h^-1))) +
  ylim(0,3)+
  ggthemes::theme_few() +
  xlab("") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

ggplot(b) +
  aes(y = interaction.index.r, x =  trt, colour = trt) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2.5, alpha = 0.7) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  ylab(bquote('Specific growth rate '(h^-1))) +
  ylim(0,3)+
  ggthemes::theme_few() +
  xlab("") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

ggplot(filter(res,Species_num>1)) +
  aes(y = interaction.r, x =  Species_num, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 1, size = 2,position=position_dodge(width=0.8)) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt),position=position_dodge(width=0.8)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt),position=position_dodge(width=0.8)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab("Interaction index") +
  xlab("Species richness") +
  ylim(0,3)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  # theme_test()+
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))


ggplot(filter(res,Species_num>1)) +
  aes(y = interaction.k, x =  Species_num, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 1, size = 2,position=position_dodge(width=0.8)) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt),position=position_dodge(width=0.8)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt),position=position_dodge(width=0.8)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab("Interaction index") +
  xlab("Species richness") +
  ylim(0,3)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  # theme_test()+
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))


# CE SE ----
res$exptd.k = c(exptd.k.ctrl,exptd.k.salt,exptd.k.stv)
res$exptd.r = c(exptd.r.ctrl,exptd.r.salt,exptd.r.stv)
res$interaction.r = res$r/res$exptd.r

res$RY.k = res$k-res$exptd.k
res$RY.r = res$r-res$exptd.r
res %>% mutate(CE.k = RY.k*exptd.k,
                    SE.k = RY.k - CE.k,
                    CE.r = RY.r*exptd.r,
                    SE.r = RY.r - CE.r) %>% filter(Species_num > 1) -> res

my_comparisons <- list(c("Control","Salinity"), c("Salinity", "Starvation"),
                       c("Control", "Starvation"))

ggplot(res) +
  aes(y = CE.k, x =  trt, colour = trt) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2, alpha = 0.8) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  ylab("CE for community yield") +
  ylim(-1.5,1.5)+
  ggthemes::theme_few() +
  xlab("") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

ggplot(res) +
  aes(y = SE.k, x =  trt, colour = trt) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2, alpha = 0.8) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  ylab(bquote('Specific growth rate '(h^-1))) +
  ylim(-0.5,0.5)+
  ggthemes::theme_few() +
  xlab("") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

# scatter 
ggplot(res) +
  aes(y = CE.k, x =  Species_num, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('Specific growth rate '(h^-1))) +
  xlab("Species richness") +
  ylim(-1.1,1)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  # theme_test()+
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))


ggplot(res) +
  aes(y = SE.k, x =  Species_num, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('Specific growth rate '(h^-1))) +
  xlab("Species richness") +
  ylim(-0.2,0.2)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  # theme_test()+
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

ggplot(res) +
  aes(y = CE.r, x =  Species_num, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('Specific growth rate '(h^-1))) +
  xlab("Species richness") +
  ylim(-2.1,2)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  # theme_test()+
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

ggplot(res) +
  aes(y = SE.r, x =  Species_num, fill = trt, colour = trt) +
  geom_point(shape = "circle", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2) + x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('Specific growth rate '(h^-1))) +
  xlab("Species richness") +
  ylim(-0.5,1)+
  ggthemes::theme_few() +
  facet_wrap(vars(trt), ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(1,2,4,8,12,16)) +
  # theme_test()+
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

