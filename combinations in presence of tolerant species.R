# combinations in presence of tolerant species
rm(list = ls())
source("code/0_data_setup.R")
comm = read.csv("../raw_data/occurrence.csv",row.names = 1) %>% filter(X1_Aeromonas_hydrophila_4AK4 + 
                                                                        X3_Escherichia_coli_DH5alpha +
                                                                        X4_Escherichia_coli_MG1655 + 
                                                                        X9_Aliiglaciecola_lipolytica_E3 +
                                                                        X13_Klebsiella_pneumoniae + 
                                                                        X14_Vibrio_cholerae +
                                                                        X15_Acinetobacter_baumannii + 
                                                                        X16_Salmonella_enterica_subsp._enterica_serovar_Typhimurium == 0) -> comm.all.tor # only 16 combinations

## combinations dominated by salt tolerance taxa
comm = read.csv("../raw_data/occurrence.csv",row.names = 1) %>% filter((X2_Cupriavidus_necator_H16 + 
                                                                         X5_Stutzerimonas_stutzeri + 
                                                                         X6_Pseudomonas_putida_KT2442 +
                                                                         X7_Proteus_vulgaris + 
                                                                         X8_Undibacterium_terreum_C3 +
                                                                         X10_Asticcacaulis_taihuensis_T3.B7 + 
                                                                         X11_Novosphingobium_taihuense_T3.B9 +
                                                                         X12_Salmonella_enterica_subsp._enterica_serovar_Braenderup)/(X1_Aeromonas_hydrophila_4AK4 + 
                                                                         X3_Escherichia_coli_DH5alpha +
                                                                         X4_Escherichia_coli_MG1655 + 
                                                                         X9_Aliiglaciecola_lipolytica_E3 +
                                                                         X13_Klebsiella_pneumoniae + 
                                                                         X14_Vibrio_cholerae +
                                                                         X15_Acinetobacter_baumannii + 
                                                                         X16_Salmonella_enterica_subsp._enterica_serovar_Typhimurium) > 1) -> comm.tor.dmn # only 16 combinations
names(comm)


res.long <- reshape2::melt(res,
                           measure.vars = c("Species_num","PD","MPD"),
                           value.name = "diversity",
                           variable.name = "diversity_type")
res.long <- res.long[which(res.long$Species_num_comb %in% row.names(comm.tor.dmn)),] %>% filter(trt == "Salinity")

ggplot(res.long) +
  aes(y = K, x =  diversity) + 
  geom_point(shape = "circle", alpha = 0.8, size = 2, color = "#6aa84f") +
  # geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
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
ggsave("plot/maximal_yield/K_diversity_tolerant_dominance.pdf", width = 7.5, height = 2.7)

fm1 = lm(K ~ diversity, filter(res.long, diversity_type == "Species_num"))
summary(fm1)
extractAIC(fm1)
fm1 = lm(K ~ diversity, filter(res.long, diversity_type == "PD"))
summary(fm1)
extractAIC(fm1)

## combinations not dominated by salt tolerance taxa
comm = read.csv("../raw_data/occurrence.csv",row.names = 1) %>% filter((X2_Cupriavidus_necator_H16 + 
                                                                          X5_Stutzerimonas_stutzeri + 
                                                                          X6_Pseudomonas_putida_KT2442 +
                                                                          X7_Proteus_vulgaris + 
                                                                          X8_Undibacterium_terreum_C3 +
                                                                          X10_Asticcacaulis_taihuensis_T3.B7 + 
                                                                          X11_Novosphingobium_taihuense_T3.B9 +
                                                                          X12_Salmonella_enterica_subsp._enterica_serovar_Braenderup)/(X1_Aeromonas_hydrophila_4AK4 + 
                                                                                                                                         X3_Escherichia_coli_DH5alpha +
                                                                                                                                         X4_Escherichia_coli_MG1655 + 
                                                                                                                                         X9_Aliiglaciecola_lipolytica_E3 +
                                                                                                                                         X13_Klebsiella_pneumoniae + 
                                                                                                                                         X14_Vibrio_cholerae +
                                                                                                                                         X15_Acinetobacter_baumannii + 
                                                                                                                                         X16_Salmonella_enterica_subsp._enterica_serovar_Typhimurium) <= 1) -> comm.sen.dmn # only 16 combinations
names(comm)


res.long <- reshape2::melt(res,
                           measure.vars = c("Species_num","PD","MPD"),
                           value.name = "diversity",
                           variable.name = "diversity_type")
res.long <- res.long[which(res.long$Species_num_comb %in% row.names(comm.sen.dmn)),] %>% filter(trt == "Salinity")

ggplot(res.long) +
  aes(y = K, x =  diversity) + 
  geom_point(shape = "circle", alpha = 0.8, size = 2, color = "#6aa84f") +
  # geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
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
ggsave("plot/maximal_yield/K_diversity_tolerant_dominance.pdf", width = 7.5, height = 2.7)

fm1 = lm(K ~ diversity, filter(res.long, diversity_type == "Species_num"))
summary(fm1)
extractAIC(fm1)
fm1 = lm(K ~ diversity, filter(res.long, diversity_type == "PD"))
summary(fm1)
extractAIC(fm1)

## combinations dominated by hungry tolerance taxa
comm = read.csv("../raw_data/occurrence.csv",row.names = 1) %>% filter((X2_Cupriavidus_necator_H16 + 
                                                                        X10_Asticcacaulis_taihuensis_T3.B7)/(
                                                                        X2_Cupriavidus_necator_H16 + 
                                                                         X5_Stutzerimonas_stutzeri + 
                                                                         X6_Pseudomonas_putida_KT2442 +
                                                                         X7_Proteus_vulgaris + 
                                                                         X8_Undibacterium_terreum_C3 +
                                                                         X10_Asticcacaulis_taihuensis_T3.B7 + 
                                                                         X11_Novosphingobium_taihuense_T3.B9 +
                                                                         X12_Salmonella_enterica_subsp._enterica_serovar_Braenderup +
                                                                           X1_Aeromonas_hydrophila_4AK4 + 
                                                                           X3_Escherichia_coli_DH5alpha +
                                                                           X4_Escherichia_coli_MG1655 +
                                                                           X9_Aliiglaciecola_lipolytica_E3 +
                                                                           X13_Klebsiella_pneumoniae + 
                                                                           X14_Vibrio_cholerae +
                                                                           X15_Acinetobacter_baumannii + 
                                                                           X16_Salmonella_enterica_subsp._enterica_serovar_Typhimurium) >= 1) -> comm.tor.dmn # only 16 combinations
names(comm)


res.long <- reshape2::melt(res,
                           measure.vars = c("Species_num","PD","MPD"),
                           value.name = "diversity",
                           variable.name = "diversity_type")
res.long <- res.long[which(res.long$Species_num_comb %in% row.names(comm.tor.dmn)),] %>% filter(trt == "Starvation")

ggplot(res.long) +
  aes(y = K, x =  diversity) + 
  geom_point(shape = "circle", alpha = 0.8, size = 2, color = "#6aa84f") +
  # geom_smooth(method = "lm", formula = y ~ x, size = 1, se = F, aes(color = trt)) +
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
ggsave("plot/maximal_yield/K_diversity_tolerant_dominance.pdf", width = 7.5, height = 2.7)

fm1 = lm(K ~ diversity, filter(res.long, diversity_type == "Species_num"))
summary(fm1)
extractAIC(fm1)
fm1 = lm(K ~ diversity, filter(res.long, diversity_type == "PD"))
summary(fm1)
extractAIC(fm1)


