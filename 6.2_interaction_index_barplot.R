# interaction index
# 为了计算richness = 16时的
rm(list = ls())
source("code/0_data_setup.R")
ccodes = c("#000000", "#6aa84f","#61a0d9")

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

comm = read.csv("../raw_data/occurrence.csv",row.names = 1)
comm$Species_num_comb = rownames(comm)

comm.ctrl = lapply(res_normal$Species_num_comb, function(x){
  row = comm[which(rownames(comm) == x),-17]
  row
})
comm.ctrl = purrr::reduce(comm.ctrl, rbind.data.frame)

comm.salt = lapply(res_salt$Species_num_comb, function(x){
  row = comm[which(rownames(comm) == x),-17]
  row
})
comm.salt = purrr::reduce(comm.salt, rbind.data.frame)

comm.stv = lapply(res_hungry$Species_num_comb, function(x){
  row = comm[which(rownames(comm) == x),-17]
  row
})
comm.stv = purrr::reduce(comm.stv, rbind.data.frame)

# data preparation ----
## growth rate ----
filter(res, Species_num == 1&trt == "Control") %>% dplyr::select(.,log_theoretical_r) -> mono.ctrl
exptd.r.ctrl = apply(comm.ctrl,1, function(x){
  sum = sum(x*mono.ctrl)
  count = sum(x)
  exp.r = sum/count
  exp.r
})

filter(res, Species_num == 1&trt == "Salinity") %>% dplyr::select(.,log_theoretical_r) -> mono.salt 
exptd.r.salt = apply(comm.salt,1, function(x){
  sum = sum(x*mono.salt)
  count = sum(x)
  exp.r = sum/count
  exp.r
})

filter(res, Species_num == 1&trt == "Starvation") %>% dplyr::select(.,log_theoretical_r) -> mono.stv
exptd.r.stv = apply(comm.stv,1, function(x){
  sum = sum(x*mono.stv)
  count = sum(x)
  exp.r = sum/count
  exp.r
})

res_full$exptd.r = c(exptd.r.ctrl,exptd.r.salt,exptd.r.stv)
res_full$interaction.r = res_full$log_theoretical_r/res_full$exptd.r

# subset species richness > 1
res_full %>% filter(Species_num >1) %>% group_by(trt, Species_num) %>% 
  summarise(interaction.index.r = mean(interaction.r), ii.sd.r = sd(interaction.r)) %>% as.data.frame() -> r.interact

## yield ----
filter(res, Species_num == 1&trt == "Control") %>% dplyr::select(.,K) -> mono.ctrl
exptd.k.ctrl = apply(comm.ctrl,1, function(x){
  sum = sum(x*mono.ctrl)
  count = sum(x)
  exp.k = sum/count
  exp.k
})

filter(res, Species_num == 1&trt == "Salinity") %>% dplyr::select(.,K) -> mono.salt 
exptd.k.salt = apply(comm.salt,1, function(x){
  sum = sum(x*mono.salt)
  count = sum(x)
  exp.k = sum/count
  exp.k
})

filter(res, Species_num == 1&trt == "Starvation") %>% dplyr::select(.,K) -> mono.stv
exptd.k.stv = apply(comm.stv,1, function(x){
  sum = sum(x*mono.stv)
  count = sum(x)
  exp.k = sum/count
  exp.k
})

res_full$exptd.k = c(exptd.k.ctrl,exptd.k.salt,exptd.k.stv)
res_full$interaction.k = res_full$K/res_full$exptd.k

# subset species richness > 1
res_full %>% filter(Species_num >1) %>% group_by(trt, Species_num) %>% 
  summarise(interaction.index.k = mean(interaction.k), ii.sd.k = sd(interaction.k)) %>% as.data.frame() -> k.interact

## death rate ----
filter(res, Species_num == 1&trt == "Control") %>% dplyr::select(.,decrease_rate) -> mono.ctrl
exptd.d.ctrl = apply(comm.ctrl,1, function(x){
  sum = sum(x*mono.ctrl)
  count = sum(x)
  exp.d = sum/count
  exp.d
})

filter(res, Species_num == 1&trt == "Salinity") %>% dplyr::select(.,decrease_rate) -> mono.salt 
exptd.d.salt = apply(comm.salt,1, function(x){
  sum = sum(x*mono.salt)
  count = sum(x)
  exp.d = sum/count
  exp.d
})

filter(res, Species_num == 1&trt == "Starvation") %>% dplyr::select(.,decrease_rate) -> mono.stv
exptd.d.stv = apply(comm.stv,1, function(x){
  sum = sum(x*mono.stv)
  count = sum(x)
  exp.d = sum/count
  exp.d
})

res_full$exptd.d = c(exptd.d.ctrl,exptd.d.salt,exptd.d.stv)
res_full$interaction.d = res_full$decrease_rate/res_full$exptd.d

# subset species richness > 1
res_full %>% filter(Species_num >1) %>% group_by(trt, Species_num) %>% 
  summarise(interaction.index.d = mean(interaction.d), ii.sd.d = sd(interaction.d)) %>% as.data.frame() -> d.interact


res_full %>% filter(Species_num >1) %>% group_by(trt, Species_num) %>% 
  dplyr::summarise(ii.r = mean(interaction.r), ii.r.sd = std.error(interaction.r),
                   ii.r.lower = ii.r - 1.96*ii.r.sd, ii.r.upper = ii.r + 1.96*ii.r.sd,
                   ii.k = mean(interaction.k), ii.k.sd = std.error(interaction.k),
                   ii.k.lower = ii.k - 1.96*ii.k.sd, ii.k.upper = ii.k + 1.96*ii.k.sd,
                   ii.d = mean(interaction.d), ii.d.sd = std.error(interaction.d),
                   ii.d.lower = ii.d - 1.96*ii.d.sd, ii.d.upper = ii.d + 1.96*ii.d.sd) %>% 
  as.data.frame() -> bar
bar$Species_num = as.factor(bar$Species_num)

# plot ----
## growth rate ----
# new plot (idea taken from Clara Qin)
for (i in 1:nrow(bar)){
  bar$ii.r.sig[i] = (bar$ii.r.lower[i] > 1 | bar$ii.r.upper[i] <1)
}
ggplot(bar, aes(x=Species_num, y=ii.r, col=trt, shape = ii.r.sig)) +#, shape=sig, alpha=sig
  geom_abline(intercept=1, slope=0, col="black", linetype = "dashed") +
  scale_shape_manual(values=c(1,16)) +
  geom_errorbar(aes(ymin=ii.r.lower, ymax=ii.r.upper), 
                position=position_dodge(0.6),
                width=0) +
  ylab("Interaction index") +
  xlab("Species richness") +
  ylim(0.75,1.85) +
  geom_point(position=position_dodge(0.6), size=2.5) + 
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
ggsave("plot/interaction/r_with_sig_barplot.pdf", width = 4.3, height = 2.9)

t1 <- aov(interaction.r ~ trt, data = filter(res_full, Species_num > 1))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.k ~ trt, data = filter(res_full, Species_num > 1))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.d ~ trt, data = filter(res_full, Species_num > 1))
duncan.test(t1,"trt",group = T, console = T)

# duncan test
library(agricolae)
t1 <- aov(interaction.r ~ trt, data = filter(res_full, Species_num ==2))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.r ~ trt, data = filter(res_full, Species_num ==4))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.r ~ trt, data = filter(res_full, Species_num ==8))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.r ~ trt, data = filter(res_full, Species_num ==12))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.r ~ trt, data = filter(res_full, Species_num ==16))
duncan.test(t1,"trt",group = T, console = T)

## yield ----
for (i in 1:nrow(bar)){
  bar$ii.k.sig[i] = (bar$ii.k.lower[i] > 1 | bar$ii.k.upper[i] <1)
}
ggplot(bar, aes(x=Species_num, y=ii.k, col=trt, shape = ii.k.sig)) +#, shape=sig, alpha=sig
  geom_abline(intercept=1, slope=0, col="black", linetype = "dashed") +
  scale_shape_manual(values=c(1,16)) +
  geom_errorbar(aes(ymin=ii.k.lower, ymax=ii.k.upper), 
                position=position_dodge(0.6),
                width=0) +
  ylab("Interaction index") +
  xlab("Species richness") +
  geom_point(position=position_dodge(0.6), size=2.5) + 
  scale_colour_manual(values=ccodes) +
  ylim(0.75,1.56) +
  ggthemes::theme_few() +
  theme(axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(color = "black", size = 12, family = "Arial"),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, family = "Arial", margin = margin(t = 8, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 0.75),
        legend.key.size = unit(0.75, 'lines'),
        legend.text = element_text(size=12, family = "Arial"),
        legend.spacing.y = unit(0.75, 'lines'))
ggsave("plot/interaction/k_with_sig_barplot.pdf", width = 4.3, height = 2.9)

# duncan test
library(agricolae)
t1 <- aov(interaction.k ~ trt, data = filter(res_full, Species_num ==2))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.k ~ trt, data = filter(res_full, Species_num ==4))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.k ~ trt, data = filter(res_full, Species_num ==8))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.k ~ trt, data = filter(res_full, Species_num ==12))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.k ~ trt, data = filter(res_full, Species_num ==16))
duncan.test(t1,"trt",group = T, console = T)

## death rate ----
for (i in 1:nrow(bar)){
  bar$ii.d.sig[i] = (bar$ii.d.lower[i] > 1 | bar$ii.d.upper[i] <1)
}
ggplot(bar, aes(x=Species_num, y=ii.d, col=trt, shape = ii.d.sig)) +#, shape=sig, alpha=sig
  geom_abline(intercept=1, slope=0, col="black", linetype = "dashed") +
  scale_shape_manual(values=c(1,16)) +
  geom_errorbar(aes(ymin=ii.d.lower, ymax=ii.d.upper), 
                position=position_dodge(0.6),
                width=0) +
  ylab("Interaction index") +
  xlab("Species richness") +
  geom_point(position=position_dodge(0.6), size=2.5) + 
  scale_colour_manual(values=ccodes) +
  ylim(0,3.5) +
  ggthemes::theme_few() +
  theme(axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(color = "black", size = 12, family = "Arial"),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", size = 12, family = "Arial", margin = margin(t = 8, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 0.75),
        legend.key.size = unit(0.75, 'lines'),
        legend.text = element_text(size=12, family = "Arial"),
        legend.spacing.y = unit(0.75, 'lines'))
ggsave("plot/interaction/d_with_sig_barplot.pdf", width = 4.3, height = 2.9)

# duncan test
library(agricolae)
t1 <- aov(interaction.d ~ trt, data = filter(res_full, Species_num ==2))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.d ~ trt, data = filter(res_full, Species_num ==4))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.d ~ trt, data = filter(res_full, Species_num ==8))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.d ~ trt, data = filter(res_full, Species_num ==12))
duncan.test(t1,"trt",group = T, console = T)
t1 <- aov(interaction.d ~ trt, data = filter(res_full, Species_num ==16))
duncan.test(t1,"trt",group = T, console = T)

# Miscellaneous ----
cor.test(res_full$k, res_full$interaction.k)
cor.test(res_full$k, res_full$copy_num)
cor.test(res_full$k, res_full$copy_num.sd)

summary(lm(interaction.k ~ Species_num*trt*MPD*copy_num*copy_num.sd, res_full))

# pca ---
res.pca = scale(dplyr::select(res_full, Species_num, MPD, copy_num, copy_num.sd, PD) %>% filter(., Species_num >1))

library(factoextra)
pca = FactoMineR::PCA(res.pca, graph = FALSE)
eig.val <- get_eigenvalue(pca)
eig.val
var <- get_pca_var(pca)
head(var$coord)
fviz_pca_var(pca, col.var = "black")
var$cos2

res.pca2 = na.omit(res.pca)
cor(res.pca2)

Dim1 = as.matrix(scale(dplyr::select(res_full, Species_num, MPD, copy_num, copy_num.sd, PD))) %*% as.matrix(var$coord[,1])
Dim2 = as.matrix(scale(dplyr::select(res_full, Species_num, MPD, copy_num, copy_num.sd, PD))) %*% as.matrix(var$coord[,2])

summary(lm(res_full$interaction.k ~ Dim1*Dim2))

summary(lm(interaction.k ~ trt*MPD*copy_num*copy_num.sd, res_full))

summary(lm(k ~ trt*Species_num, res_full))
summary(lm(r ~ trt*Species_num, res_full))

summary(lm(k ~ trt*Species_num*exptd.k, filter(res_full, Species_num > 1)))
anova(lm(k ~ trt*Species_num*exptd.k, filter(res_full, Species_num > 1)))
summary(lm(r ~ trt*Species_num*exptd.r, filter(res_full, Species_num > 1)))

summary(lm(interaction.r ~ trt*Species_num, filter(res_full, Species_num > 1)))

summary(lm(r ~ Species_num*exptd.r, filter(res_full, Species_num > 1,trt == "Starvation")))
summary(lm(r ~ Species_num*exptd.r, filter(res_full, Species_num > 1,trt == "Salinity")))
summary(lm(r ~ Species_num*exptd.r, filter(res_full,Species_num > 1, trt == "Control")))


summary(lm(k ~ Species_num*exptd.k*trt, filter(res_full, Species_num > 1)))
summary(lm(r ~ Species_num*exptd.r, filter(res_full, Species_num > 1,trt == "Salinity")))
summary(lm(r ~ Species_num*exptd.r, filter(res_full,Species_num > 1, trt == "Control")))

anova(lm(r ~ trt*Species_num*exptd.r, res_full))


cor.test(filter(res_full, Species_num > 1, trt == "Salinity")$interaction.r,  filter(res_full, Species_num > 1, trt == "Salinity")$Species_num)
cor.test(filter(res_full, Species_num > 1, trt == "Starvation")$interaction.r,  filter(res_full, Species_num > 1, trt == "Starvation")$Species_num)

anova(lm(interaction.k ~ trt*MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1)))
anova(lm(interaction.k ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Control")))
anova(lm(interaction.k ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Salinity")))
anova(lm(interaction.k ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Starvation")))

summary(lm(interaction.k ~ trt*MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1)))

anova(lm(interaction.r ~ trt*MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1)))
anova(lm(interaction.r ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Control")))
anova(lm(interaction.r ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Salinity")))
anova(lm(interaction.r ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Starvation")))

summary(lm(interaction.r ~ trt*MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1)))

summary(lm(r ~ trt*MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1)))
summary(lm(r ~ trt*MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1)))
anova(lm(k ~ trt*MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1)))



cor.test(filter(res_full, Species_num > 1)$copy_num, filter(res_full, Species_num > 1)$interaction.r)
cor.test(filter(res_full, Species_num > 1)$MPD, filter(res_full, Species_num > 1)$interaction.k)
cor.test(filter(res_full, Species_num > 1)$copy_num.sd, filter(res_full, Species_num > 1)$interaction.k)

boxplot(res_full$r ~ res_full$trt)
boxplot(res_full$k ~ res_full$trt)

summary(lm(r ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Control")))
summary(lm(r ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Salinity")))
summary(lm(r ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Starvation")))

summary(lm(k ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Control")))
summary(lm(k ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Salinity")))
summary(lm(k ~ MPD*copy_num*copy_num.sd, filter(res_full, Species_num > 1, trt == "Starvation")))

# 不对劲 走歪了
library(agricolae)
t1 <- aov(interaction.r ~ Species_num, data = filter(res_full, Species_num > 1, trt == "Control"))
duncan.test(t1,"Species_num",group = T, console = T)

t1 <- aov(interaction.r ~ Species_num, data = filter(res_full, Species_num > 1, trt == "Salinity"))
duncan.test(t1,"Species_num",group = T, console = T)

t1 <- aov(interaction.r ~ Species_num, data = filter(res_full, Species_num > 1, trt == "Starvation"))
duncan.test(t1,"Species_num",group = T, console = T)


# ----
comm = read.csv("raw_data/occurrence.csv")

res_full = merge(res_full, comm, by = "Species_num_comb")
full.model <- lm(r ~ trt*Species_num*X1_Aeromonas_hydrophila_4AK4 +  trt*Species_num*X2_Ralstonia_eutropha_H16 + trt*Species_num*X3_Escherichia_coli_DH5 + trt*Species_num*X4_Escherichia_coli_MG1655 +
                  trt*Species_num*X5_Pseudomonas_stutzeri_1317 + trt*Species_num*X6_Pseudomonas_putida_KT2440 +
                  trt*Species_num*X7_Proteus_vulgaris +  trt*Species_num*X8_Undibacterium_terreum+trt*Species_num*X9_Glaciecola_lipolytica +
                  trt*Species_num*X10_Asticcacaulis_taihuensis + trt*Species_num*X11_Novosphingobium_taihuense +   
                  trt*Species_num*X12_Salmonella_Braenderup_serotype_H9812 +  trt*Species_num*X13_Klebsiella_pneumoniae + trt*Species_num*X14_Vibrio_cholera +
                   trt*Species_num*X15_Acinetobacter_baumannii +  trt*Species_num*X16_Salmonella_Typhimurium_ATCC_14028, data = res_full)
step.model <- stepAIC(full.model, direction = "both", 
                      trace = T)
summary(step.model)
anova(step.model)
car::Anova(lm, type = 2, singular.ok = T)

full.model <- lm(k ~ Species_num*X1_Aeromonas_hydrophila_4AK4 +  Species_num*X2_Ralstonia_eutropha_H16 + Species_num*X3_Escherichia_coli_DH5 + Species_num*X4_Escherichia_coli_MG1655 +
                  Species_num*X5_Pseudomonas_stutzeri_1317 + Species_num*X6_Pseudomonas_putida_KT2440 +
                   Species_num*X7_Proteus_vulgaris + Species_num*X8_Undibacterium_terreum+Species_num*X9_Glaciecola_lipolytica +
                   Species_num*X10_Asticcacaulis_taihuensis + Species_num*X11_Novosphingobium_taihuense +   
                   Species_num*X12_Salmonella_Braenderup_serotype_H9812 +  Species_num*X13_Klebsiella_pneumoniae + Species_num*X14_Vibrio_cholera +
                   Species_num*X15_Acinetobacter_baumannii +  Species_num*X16_Salmonella_Typhimurium_ATCC_14028, data = filter(res_full, trt == "Control"))
step.model <- stepAIC(full.model, direction = "both", 
                      trace = T)
summary(step.model)
anova(step.model)

full.model <- lm(k ~ Species_num*X1_Aeromonas_hydrophila_4AK4 +  Species_num*X2_Ralstonia_eutropha_H16 + Species_num*X3_Escherichia_coli_DH5 + Species_num*X4_Escherichia_coli_MG1655 +
                   Species_num*X5_Pseudomonas_stutzeri_1317 + Species_num*X6_Pseudomonas_putida_KT2440 +
                   Species_num*X7_Proteus_vulgaris + Species_num*X8_Undibacterium_terreum+Species_num*X9_Glaciecola_lipolytica +
                   Species_num*X10_Asticcacaulis_taihuensis + Species_num*X11_Novosphingobium_taihuense +   
                   Species_num*X12_Salmonella_Braenderup_serotype_H9812 +  Species_num*X13_Klebsiella_pneumoniae + Species_num*X14_Vibrio_cholera +
                   Species_num*X15_Acinetobacter_baumannii +  Species_num*X16_Salmonella_Typhimurium_ATCC_14028, data = filter(res_full, trt == "Salinity"))
step.model <- stepAIC(full.model, direction = "both", 
                      trace = T)
summary(step.model)
anova(step.model)

full.model <- lm(k ~ Species_num*X1_Aeromonas_hydrophila_4AK4 +  Species_num*X2_Ralstonia_eutropha_H16 + Species_num*X3_Escherichia_coli_DH5 + Species_num*X4_Escherichia_coli_MG1655 +
                   Species_num*X5_Pseudomonas_stutzeri_1317 + Species_num*X6_Pseudomonas_putida_KT2440 +
                   Species_num*X7_Proteus_vulgaris + Species_num*X8_Undibacterium_terreum+Species_num*X9_Glaciecola_lipolytica +
                   Species_num*X10_Asticcacaulis_taihuensis + Species_num*X11_Novosphingobium_taihuense +   
                   Species_num*X12_Salmonella_Braenderup_serotype_H9812 +  Species_num*X13_Klebsiella_pneumoniae + Species_num*X14_Vibrio_cholera +
                   Species_num*X15_Acinetobacter_baumannii +  Species_num*X16_Salmonella_Typhimurium_ATCC_14028, data = filter(res_full, trt == "Starvation"))
step.model <- stepAIC(full.model, direction = "both", 
                      trace = T)
summary(step.model)
anova(step.model)


c = res_full %>% filter(.,Species_num >1, trt == "Control") %>% as.data.frame(.)
t.test(c$interaction.r, mu = 1, alternative = "two.sided")
std.error(c$interaction.r)

c = res_full %>% filter(.,Species_num >1, trt == "Salinity") %>% as.data.frame(.)
t.test(c$interaction.r, mu = 1, alternative = "two.sided")
std.error(c$interaction.r)

c = res_full %>% filter(.,Species_num >1, trt == "Starvation") %>% as.data.frame(.)
t.test(c$interaction.r, mu = 1, alternative = "two.sided")
std.error(c$interaction.r)


c = res_full %>% filter(.,Species_num >1, trt == "Control") %>% as.data.frame(.)
t.test(c$interaction.k, mu = 1, alternative = "two.sided")
std.error(c$interaction.k)

c = res_full %>% filter(.,Species_num >1, trt == "Salinity") %>% as.data.frame(.)
t.test(c$interaction.k, mu = 1, alternative = "two.sided")
std.error(c$interaction.k)

c = res_full %>% filter(.,Species_num >1, trt == "Starvation") %>% as.data.frame(.)
t.test(c$interaction.k, mu = 1, alternative = "two.sided")
std.error(c$interaction.k)


