# rrn comparison
source("code/0_data_setup.R")
rrn =  read.csv("data/rrn_tolerance/rrn.csv",row.names = 1)

mean(rrn$rrn.copy.number) # 5.75
sum(rrn$rrn.copy.number > 5.75)
rrn$rrn.copy.number > 5.75 # 1 4 6 7 12 13 14 15 16

comm = read.csv("../raw_data/occurrence.csv",row.names = 1) 
names(comm)

comm[1:16,] %>% filter(X1_Aeromonas_hydrophila_4AK4 + 
             X3_Escherichia_coli_DH5alpha +
             X4_Escherichia_coli_MG1655 + 
             X9_Aliiglaciecola_lipolytica_E3 +
             X13_Klebsiella_pneumoniae + 
             X14_Vibrio_cholerae +
             X15_Acinetobacter_baumannii + 
             X16_Salmonella_enterica_subsp._enterica_serovar_Typhimurium == 0) -> comm.all.tor

res.annotated <- filter(res, Species_num == 1) %>% mutate(rrn.level = ifelse(copy_num > 5.75,"high","low")) 

ggplot(res.annotated) +
  aes(y = log_theoretical_r, x =  rrn.level) +
  # geom_boxplot(outlier.shape = NA) +
  geom_violin(aes(fill = rrn.level), trim = FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 1, alpha = 0.5) + 
  geom_boxplot(width = 0.15) +
  # geom_boxplot(width=.1)+
  stat_compare_means(method = "t.test", paired = F, label.y = 0.04 ) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  ylab(bquote('Maximal specific growth rate '(min^-1))) +
  # ylim(0,0.045)+
  ggthemes::theme_few() +
  xlab("") +
  facet_wrap("trt")+
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

ggplot(res.annotated) +
  aes(y = K, x =  rrn.level) +
  # geom_boxplot(outlier.shape = NA) +
  geom_violin(aes(fill = rrn.level), trim = FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 1, alpha = 0.5) + 
  geom_boxplot(width = 0.15) +
  # geom_boxplot(width=.1)+
  stat_compare_means(method = "t.test", paired = F, label.y = 2.5) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  ylab(bquote('Maximal yield (OD units)')) +
  #ylim(0,0.045)+
  ggthemes::theme_few() +
  xlab("") +
  facet_wrap("trt")+
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

ggplot(res.annotated) +
  aes(y = decrease_rate, x =  rrn.level) +
  # geom_boxplot(outlier.shape = NA) +
  geom_violin(aes(fill = rrn.level), trim = FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 1, alpha = 0.5) + 
  geom_boxplot(width = 0.15) +
  # geom_boxplot(width=.1)+
  stat_compare_means(method = "t.test", paired = F,label.y = 7) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  ylab(bquote('Death rate '(~x10^-4 ~ min^-1))) +
  #ylim(0,0.045)+
  ggthemes::theme_few() +
  xlab("") +
  facet_wrap("trt")+
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))

# all combinations
ggplot(res) +
  aes(y = decrease_rate, x =  copy_num, fill = trt, colour = trt) + 
  geom_point(shape = "circle", alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
  geom_smooth(method = "lm", formula = y ~ I(x^2)+x, size = 1, se = F, aes(color = trt)) +
  scale_fill_manual(values = ccodes) +
  scale_color_manual(values = ccodes) +
  expand_limits(x = 0) +
  ylab(bquote('Death rate '(~x10^-4 ~ min^-1))) +
  xlab("Biodiversity") +
  # ylim(0,10)+
  # ggthemes::theme_few() +
  facet_wrap("trt", scales = "free", strip.position = "top") +
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
#              parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5) +
# stat_poly_eq(formula = y ~ I(x^2)+x, 
#              aes(label = paste(..eq.label..,..p.value.label.., sep = "~~~~")), 
#              parse = TRUE,coef.digits = 3,p.digits = 3, coef.keep.zeros = F, size = 2.5,
#              vjust = 1.5)