# Linear models followed by ANOVA
rm(list = ls())
source("code/0_data_setup.R")
comm = read.csv("../raw_data/occurrence.csv")
df = merge(res, comm, by = "Species_num_comb")

# Interaction of trt and diversity ----
sink("data/linear_models/trt_diversity.txt")
print("---- Species richness ----")
lm <- lm(log_theoretical_r ~ Species_num*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

lm <- lm(K ~ Species_num*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

lm <- lm(decrease_rate ~ Species_num*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

print("---- Faith's PD ----")
lm <- lm(log_theoretical_r ~ PD*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

lm <- lm(K ~ PD*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

lm <- lm(decrease_rate ~ PD*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

print("---- MPD ----")
lm <- lm(log_theoretical_r ~ MPD*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

lm <- lm(K ~ MPD*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

lm <- lm(decrease_rate ~ MPD*trt, data = df)
summary(lm)
car::Anova(lm, type = 2)

closeAllConnections()

# Interaction of trt and diversity ----
sink("data/linear_models/tol_species.txt")
print("---- Species richness ----")



