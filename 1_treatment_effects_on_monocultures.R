# treatment effects on different species in monocultures
rm(list = ls())
source("code/0_data_setup.R")

# Treatment effects on monocultures
read.csv("data/mono/mono_data_statistical_test_R.csv", row.names = 1) -> df
star_position = function(x){
  id_posi = which(x>0)
  id_neg = which(x<0)
  x[id_posi] = x[id_posi] + 7.5
  x[id_neg] = x[id_neg] - 7.5
  return(x)
}

star_p = star_position(df$perc_change_in_K)

ggplot(df, aes(y = perc_change_in_K, x = Species.ID,  group = trt, fill=trt)) +
  geom_bar(stat="identity", width = 0.75, position=position_dodge(), color = "black", alpha = 0.6, size = 0.3) +
  geom_errorbar(aes(ymin=perc_change_in_K-perc_change_in_K.se, ymax=perc_change_in_K+perc_change_in_K.se), width=.2, position=position_dodge(.8)) +
  theme_test() +
  ylab("Percentage change in maximal yield (%)") +
  # theme_test()+
  theme(strip.text = element_text(family = "Arial",size = 12),
        axis.text = element_text(color = "black", size = 12, family = "Arial"),
        axis.title.y = element_text(family = "Arial", color = "black", size = 12, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        # legend.position = "none",
        axis.title.x = element_text(color = "black", family = "Arial", size = 12, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1)) +
  # scale_fill_manual(values=c('black','lightgray')) +
  geom_hline(yintercept = 0) +
  geom_text(aes(y=star_p, label = K.sig.label),position = position_dodge(0.8),size=5)+
  xlab("Species ID") +
  scale_fill_manual(values = ccodes[2:3])
ggsave("plot/trt_effects_on_mono/K.pdf", width = 7.5, height = 4.2)

