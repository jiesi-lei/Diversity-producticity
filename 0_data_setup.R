# load packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(egg)
library(RColorBrewer)
library(MASS)
library(sfsmisc)
library(ggsci)
library("scales")
library(wesanderson)
library(plotrix)
library(ggh4x)#set scales for facet_wrap

library(showtext)
font_add("Arial", "/Library/Fonts/Arial.ttf")  # Use the actual file path
showtext_auto()

theme_set(theme_test(base_family = "Arial"))

# load data
group_normal = read.csv("../raw_data/normal_sample_list_clear.csv",row.names = 1,stringsAsFactors = T)
res_normal = read.csv("results_matlab/Feb07/result-normal-final-Feb07.csv",row.names = 1)
res_normal = res_normal[which(res_normal$Sample %in% group_normal$Original_name)[1:nrow(group_normal)],]
sum(res_normal$Sample != group_normal$Original_name) # 0 indicates that they matched well
res_normal %>% mutate(Rmax = -log_theoretical_r*log_theoretical_K/4) %>%
  mutate(Species_num = group_normal$Species_num,
         Combination = group_normal$Combination,
         Species_num_comb = group_normal$Species_num_comb,
         decrease_rate = -decrease_rate*10000) %>% 
  #filter(.,decrease_rate > 0) %>% 
  na.omit(.) -> res_normal 


group_salt = read.csv("../raw_data/salt_sample_list_clear.csv",row.names = 1,stringsAsFactors = T)
res_salt = read.csv("results_matlab/Feb07/result-salt-final-Feb07.csv",row.names = 1)
#which(res_salt$Sample != group_salt$Original_name)
res_salt = res_salt[which(res_salt$Sample %in% group_salt$Original_name)[1:nrow(group_salt)],]
sum(res_salt$Sample != group_salt$Original_name) # 0 indicates that they matched well
res_salt %>% mutate(Rmax = -log_theoretical_r*log_theoretical_K/4) %>%
  mutate(Species_num = group_salt$Species_num,
         Combination = group_salt$Combination,
         Species_num_comb = group_salt$Species_num_comb,
         decrease_rate = -decrease_rate*10000) %>% 
  #filter(.,decrease_rate < 1e-3 & decrease_rate > 0)%>% 
  na.omit(.) -> res_salt

group_hungry = read.csv("../raw_data/hungry_sample_list_clear.csv",row.names = 1,stringsAsFactors = T)
res_hungry = read.csv("results_matlab/Feb07/result-hungry-final-Feb07.csv",row.names = 1)
res_hungry = res_hungry[which(res_hungry$Sample %in% group_hungry$Original_name)[1:nrow(group_hungry)],]
#which(res_hungry$Sample != group_hungry$Original_name)
sum(res_hungry$Sample != group_hungry$Original_name) # 0 indicates that they matched well
res_hungry %>% mutate(Rmax = -log_theoretical_r*log_theoretical_K/4) %>%
  mutate(Species_num = group_hungry$Species_num,
         Combination = group_hungry$Combination,
         Species_num_comb = group_hungry$Species_num_comb,
         decrease_rate = -decrease_rate*10000)  %>% 
  #filter(.,decrease_rate > 0) %>% 
  na.omit(.) -> res_hungry

# remove outliers and average across replicates
out <- boxplot(res_normal$log_theoretical_r)$out
out_ind <- which(res_normal$log_theoretical_r %in% out)
out_ind
# 55,57,71,79,88,107,116,207
res_normal = res_normal[-c(37,89), ]
# res_normal = res_normal[-c(which(res_normal$exponential_coef_simple == 0 | res_normal$log_theoretical_K == 0)),]
aggregate(res_normal[,2:35], list(res_normal$Species_num_comb), mean) %>% mutate(across(2:35, round, 6)) -> res_normal_mean # average across replicates
aggregate(res_normal[,2:35], list(res_normal$Species_num_comb), std.error) %>% mutate(across(2:35, round, 6)) -> res_normal_se
names(res_normal_se)<-paste0(names(res_normal_se),".se")
group_normal = unique(group_normal[,3:5])
res_normal$log_theoretical_r = -res_normal$log_theoretical_r
aggregate(res_normal[,2:35], list(res_normal$Species_num_comb), mean) %>% mutate(across(2:35, round, 6)) %>% 
  dplyr::select(Group.1, Rmax, exponential_coef_simple, exponential_r, log_theoretical_r, K,mean_window_K, log_theoretical_K, decrease_rate) -> res_normal_mean # average across replicates
aggregate(res_normal[,2:35], list(res_normal$Species_num_comb), sd) %>% mutate(across(2:35, round, 6)) %>% 
  dplyr::select(Group.1, Rmax, exponential_coef_simple, exponential_r, log_theoretical_r, K, mean_window_K, log_theoretical_K,decrease_rate) -> res_normal_sd
res_normal_cv <- res_normal_sd[,2:9]/res_normal_mean[,2:9]
names(res_normal_cv)<-paste0(names(res_normal_cv),".cv")
group_normal <- group_normal[match(res_normal_mean$Group.1,group_normal$Species_num_comb),]
normal <- data.frame(group_normal,res_normal_mean,res_normal_se,res_normal_cv)
rm(group_normal,res_normal_cv,res_normal_mean,res_normal_sd)

out <- boxplot(res_salt$Rmax)$out
out_ind <- which(res_salt$Rmax %in% c(out))
out_ind
#5,8,187,188
res_salt = res_salt[-c(8,187,188), ]
# res_salt = res_salt[-c(which(res_salt$exponential_coef_simple == 0 | res_salt$log_theoretical_K == 0)),]
aggregate(res_salt[,2:35], list(res_salt$Species_num_comb), mean) %>% mutate(across(2:35, round, 6)) -> res_salt_mean # average across replicates
aggregate(res_salt[,2:35], list(res_salt$Species_num_comb), std.error) %>% mutate(across(2:35, round, 6)) -> res_salt_se
names(res_salt_se)<-paste0(names(res_salt_se),".se")
group_salt = unique(group_salt[,3:5])
res_salt$log_theoretical_r = -res_salt$log_theoretical_r
aggregate(res_salt[,2:35], list(res_salt$Species_num_comb), mean) %>% mutate(across(2:35, round, 6)) %>% 
  dplyr::select(Group.1, Rmax, exponential_coef_simple, exponential_r, log_theoretical_r, K, mean_window_K, log_theoretical_K, decrease_rate) -> res_salt_mean # average across replicates
aggregate(res_salt[,2:35], list(res_salt$Species_num_comb), sd) %>% mutate(across(2:35, round, 6)) %>% 
  dplyr::select(Group.1, Rmax, exponential_coef_simple, exponential_r, log_theoretical_r, K,mean_window_K, log_theoretical_K, decrease_rate) -> res_salt_sd
res_salt_cv <- res_salt_sd[,2:9]/res_salt_mean[,2:9]
names(res_salt_cv)<-paste0(names(res_salt_cv),".cv")
group_salt <- group_salt[match(res_salt_mean$Group.1,group_salt$Species_num_comb),]
salt <- data.frame(group_salt,res_salt_mean,res_salt_se,res_salt_cv)
rm(group_salt,res_salt_cv,res_salt_mean,res_salt_sd)

out <- boxplot(res_hungry$Rmax)$out
out_ind <- which(res_hungry$Rmax %in% out)
out_ind
# 118 128 143 217 238 313
#res_hungry = res_hungry[-c( 118,128,143,217,238,313), ]
# res_hungry = res_hungry[-c(which(res_hungry$exponential_coef_simple == 0 | res_hungry$log_theoretical_K == 0)),]
aggregate(res_hungry[,2:35], list(res_hungry$Species_num_comb), mean) %>% mutate(across(2:35, round, 6)) -> res_hungry_mean # average across replicates
aggregate(res_hungry[,2:35], list(res_hungry$Species_num_comb), std.error) %>% mutate(across(2:35, round, 6)) -> res_hungry_se
names(res_hungry_se)<-paste0(names(res_hungry_se),".se")
group_hungry = unique(group_hungry[,3:5])
res_hungry$log_theoretical_r = -res_hungry$log_theoretical_r
aggregate(res_hungry[,2:35], list(res_hungry$Species_num_comb), mean) %>% mutate(across(2:35, round, 6)) %>% 
  dplyr::select(Group.1, Rmax, exponential_coef_simple, exponential_r, log_theoretical_r, K, mean_window_K, log_theoretical_K,decrease_rate) -> res_hungry_mean # average across replicates
aggregate(res_hungry[,2:35], list(res_hungry$Species_num_comb), sd) %>% mutate(across(2:35, round, 6)) %>% 
  dplyr::select(Group.1, Rmax, exponential_coef_simple, exponential_r, log_theoretical_r, K, mean_window_K, log_theoretical_K,decrease_rate) -> res_hungry_sd
res_hungry_cv <- res_hungry_sd[,2:9]/res_hungry_mean[,2:9]
names(res_hungry_cv)<-paste0(names(res_hungry_cv),".cv")
group_hungry <- group_hungry[match(res_hungry_mean$Group.1,group_hungry$Species_num_comb),]
hungry <- data.frame(group_hungry,res_hungry_mean,res_hungry_se,res_hungry_cv)
rm(group_hungry,res_hungry_cv,res_hungry_mean,res_hungry_sd)


# assign PD to species combinations
PD_df = read.csv("results_R/pd_mpd_Jan2022.csv",row.names = 1)

normal$MPD <- sapply(1:nrow(normal),function(i){
  PD_df$mpd.dist[which(rownames(PD_df) == normal$Species_num_comb[i])]
})
salt$MPD <- sapply(1:nrow(salt),function(i){
  PD_df$mpd.dist[which(rownames(PD_df) == salt$Species_num_comb[i])]
})
hungry$MPD <- sapply(1:nrow(hungry),function(i){
  PD_df$mpd.dist[which(rownames(PD_df) == hungry$Species_num_comb[i])]
})

normal$PD <- sapply(1:nrow(normal),function(i){
  PD_df$PD[which(rownames(PD_df) == normal$Species_num_comb[i])]
})
salt$PD <- sapply(1:nrow(salt),function(i){
  PD_df$PD[which(rownames(PD_df) == salt$Species_num_comb[i])]
})
hungry$PD <- sapply(1:nrow(hungry),function(i){
  PD_df$PD[which(rownames(PD_df) == hungry$Species_num_comb[i])]
})

# assign average rrn copy number to species combinations
rrn = read.csv("data/rrn_tolerance/copy_number.csv",row.names = 1)

normal$copy_num <- sapply(1:nrow(normal),function(i){
  rrn$copy_number[which(rrn$rownames.comm. == normal$Species_num_comb[i])]
})
salt$copy_num <- sapply(1:nrow(salt),function(i){
  rrn$copy_number[which(rrn$rownames.comm. == salt$Species_num_comb[i])]
})
hungry$copy_num <- sapply(1:nrow(hungry),function(i){
  rrn$copy_number[which(rrn$rownames.comm. == hungry$Species_num_comb[i])]
})

normal$copy_num.sd <- sapply(1:nrow(normal),function(i){
  rrn$copy_number.sd[which(rrn$rownames.comm. == normal$Species_num_comb[i])]
})
salt$copy_num.sd <- sapply(1:nrow(salt),function(i){
  rrn$copy_number.sd[which(rrn$rownames.comm. == salt$Species_num_comb[i])]
})
hungry$copy_num.sd <- sapply(1:nrow(hungry),function(i){
  rrn$copy_number.sd[which(rrn$rownames.comm. == hungry$Species_num_comb[i])]
})

# assign number of tolerant species to species combinations
tol = read.csv("data/rrn_tolerance/tolerance.csv",row.names = 1)

normal$salt.tol <- sapply(1:nrow(normal),function(i){
  tol$salt.tol[which(rrn$rownames.comm. == normal$Species_num_comb[i])]
})
salt$salt.tol <- sapply(1:nrow(salt),function(i){
  tol$salt.tol[which(rrn$rownames.comm. == salt$Species_num_comb[i])]
})
hungry$salt.tol <- sapply(1:nrow(hungry),function(i){
  tol$salt.tol[which(rrn$rownames.comm. == hungry$Species_num_comb[i])]
})

normal$starve.tol <- sapply(1:nrow(normal),function(i){
  tol$starve.tol[which(rrn$rownames.comm. == normal$Species_num_comb[i])]
})
salt$starve.tol <- sapply(1:nrow(salt),function(i){
  tol$starve.tol[which(rrn$rownames.comm. == salt$Species_num_comb[i])]
})
hungry$starve.tol <- sapply(1:nrow(hungry),function(i){
  tol$starve.tol[which(rrn$rownames.comm. == hungry$Species_num_comb[i])]
})

normal %>% mutate(trt = "Control") -> normal
salt %>% mutate(trt = "Salinity") -> salt
hungry %>% mutate(trt = "Starvation") -> hungry
res = rbind(normal, salt, hungry)

# assign to original data
res_normal$MPD <- sapply(1:nrow(res_normal),function(i){
  PD_df$mpd.dist[which(rownames(PD_df) == res_normal$Species_num_comb[i])]
})
res_salt$MPD <- sapply(1:nrow(res_salt),function(i){
  PD_df$mpd.dist[which(rownames(PD_df) == res_salt$Species_num_comb[i])]
})
res_hungry$MPD <- sapply(1:nrow(res_hungry),function(i){
  PD_df$mpd.dist[which(rownames(PD_df) == res_hungry$Species_num_comb[i])]
})

res_normal$PD <- sapply(1:nrow(res_normal),function(i){
  PD_df$PD[which(rownames(PD_df) == res_normal$Species_num_comb[i])]
})
res_salt$PD <- sapply(1:nrow(res_salt),function(i){
  PD_df$PD[which(rownames(PD_df) == res_salt$Species_num_comb[i])]
})
res_hungry$PD <- sapply(1:nrow(res_hungry),function(i){
  PD_df$PD[which(rownames(PD_df) == res_hungry$Species_num_comb[i])]
})


# assign average rrn copy number to species combinations
res_normal$copy_num <- sapply(1:nrow(res_normal),function(i){
  rrn$copy_number[which(rrn$rownames.comm. == res_normal$Species_num_comb[i])]
})
res_salt$copy_num <- sapply(1:nrow(res_salt),function(i){
  rrn$copy_number[which(rrn$rownames.comm. == res_salt$Species_num_comb[i])]
})
res_hungry$copy_num <- sapply(1:nrow(res_hungry),function(i){
  rrn$copy_number[which(rrn$rownames.comm. == res_hungry$Species_num_comb[i])]
})

res_normal$copy_num.sd <- sapply(1:nrow(res_normal),function(i){
  rrn$copy_number.sd[which(rrn$rownames.comm. == res_normal$Species_num_comb[i])]
})
res_salt$copy_num.sd <- sapply(1:nrow(res_salt),function(i){
  rrn$copy_number.sd[which(rrn$rownames.comm. == res_salt$Species_num_comb[i])]
})
res_hungry$copy_num.sd <- sapply(1:nrow(res_hungry),function(i){
  rrn$copy_number.sd[which(rrn$rownames.comm. == res_hungry$Species_num_comb[i])]
})

# assign number of tolerant species to species combinations
res_normal$salt.tol <- sapply(1:nrow(res_normal),function(i){
  tol$salt.tol[which(rrn$rownames.comm. == res_normal$Species_num_comb[i])]
})
res_salt$salt.tol <- sapply(1:nrow(res_salt),function(i){
  tol$salt.tol[which(rrn$rownames.comm. == res_salt$Species_num_comb[i])]
})
res_hungry$salt.tol <- sapply(1:nrow(res_hungry),function(i){
  tol$salt.tol[which(rrn$rownames.comm. == res_hungry$Species_num_comb[i])]
})

res_normal$starve.tol <- sapply(1:nrow(res_normal),function(i){
  tol$starve.tol[which(rrn$rownames.comm. == res_normal$Species_num_comb[i])]
})
res_salt$starve.tol <- sapply(1:nrow(res_salt),function(i){
  tol$starve.tol[which(rrn$rownames.comm. == res_salt$Species_num_comb[i])]
})
res_hungry$starve.tol <- sapply(1:nrow(res_hungry),function(i){
  tol$starve.tol[which(rrn$rownames.comm. == res_hungry$Species_num_comb[i])]
})

res_normal %>% mutate(trt = "Control") -> res_normal
res_salt %>% mutate(trt = "Salinity") -> res_salt
res_hungry %>% mutate(trt = "Starvation") -> res_hungry
res_full = rbind(res_normal, res_salt, res_hungry)
rm(PD_df,rrn,tol)

makeStars <- function(x){
  stars <- c("***", "**", "*")
  vec <- c(0, 0.001, 0.01, 0.05)
  i <- findInterval(x, vec)
  stars[i]
}

ccodes = c("#808080", "#6aa84f","#61a0d9")
