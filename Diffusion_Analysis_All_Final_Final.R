#==============================================================================
# FIXED COPY (v5) of Diffusion_Analysis_All_Final.R  (every change is tagged #EDIT)
# To run: put this file and the six data files in ~/Downloads, then
#   source("~/Downloads/Diffusion_Analysis_All_Final_v5.R")
# When it finishes: the model tables print in the console, the plots show in
# the Plots pane, and all files are saved to ~/Downloads/DTI_results
# (ModelRes.xlsx, BehavCorr.xlsx and four PNGs). The folder path prints last.
#  1. Paths: paste(DIR, "file") joins with a space ("Downloads file.xlsx"), so
#     no file can be found. Now file.path(DIR, "file"); DIR is "~/Downloads".
#  2. ppcor: if the package isn't installed, pcor.test() falls back to a
#     base-R copy of ppcor 1.1's formula (identical results).
#  3. Saving: the commented-out write.xlsx()/png() calls are switched on.
#  4. Statistics fix, inside the six model functions only (H1, H1 Exploratory
#     and JHU; FA and MD). The old code summarized the model's fitted values
#     (predict on the observed data), so AdjMean was the raw group mean and
#     AdjSD reflected only age. Now:
#       AdjMean = age-adjusted mean: each group's model prediction at the
#                 sample's mean age (the estimated marginal mean)
#       AdjSD   = observed SD of the outcome in each group
#     AdjSMD and every line after the functions are unchanged.
#  5. Age-adjusted FA graph (FAp2; Fig. 3B): it plotted fitted values; it now
#     plots each observed FA value moved to the sample's mean age. Its y-axis
#     now matches the raw FA graph, because 13 corrected values fall outside
#     the old 0.36-0.44 axis. Nothing else about any graph changed.
#  6. Printing: model tables print in the console and plots display when the
#     script is run with source().
#  7. Gridlines: horizontal gridlines are drawn only at the labeled y-axis
#     values (FA: 0.28, 0.32 ... 0.52; MD: 7.2, 7.6 ... 9.2). Tick marks, labels,
#     axis ranges and everything else in the graphs are unchanged.
#==============================================================================

library(data.table)
library(tidyverse)
library(readr)
library(openxlsx)
library(readxl)
library(lubridate)
if (requireNamespace("ppcor", quietly = TRUE)) library(ppcor) else { #EDIT (2)
  pcor.test <- function(x, y, z, method = c("pearson", "kendall", "spearman")) {
    # same algorithm as ppcor::pcor.test (v1.1): partial r from inverse covariance
    method <- match.arg(method)
    x <- c(x); y <- c(y); z <- as.data.frame(z)
    xyz <- data.frame(x, y, z)
    n <- dim(xyz)[1]; gp <- dim(xyz)[2] - 2
    icvx <- solve(cov(xyz, method = method))
    pcor <- -cov2cor(icvx)[1, 2]
    if (method == "kendall") {
      statistic <- pcor / sqrt(2 * (2 * (n - gp) + 5) / (9 * (n - gp) * (n - 1 - gp)))
      p.value <- 2 * pnorm(-abs(statistic))
    } else {
      statistic <- pcor * sqrt((n - 2 - gp) / (1 - pcor^2))
      p.value <- 2 * pt(-abs(statistic), (n - 2 - gp))
    }
    data.frame(estimate = pcor, p.value = p.value, statistic = statistic,
               n = n, gp = gp, Method = method)
  }
}
library(ggplot2)
library(ggpubr)
library(car)

DIR<-"~/Downloads" #EDIT (1)
OUT<-file.path(DIR, "DTI_results"); dir.create(OUT, showWarnings = FALSE) #EDIT (3)
#Jasaon DTI paper data
data<-read_excel(file.path(DIR, "TensorBO_wm_FA_MD_BAS_km.xlsx"), na = "NA") #EDIT (1)
data <- data %>%
  mutate(across(3:20, as.numeric))

data2<-read_excel(file.path(DIR, "Combined_master101624_redcap_v1.xlsx"),  #EDIT (1)
                  sheet = "Combined_master101624_clean")
demo<-data2[,c("Subject ID", "pain_years_num", "age")];colnames(demo)[1]<-"subject"
demo$subject<-gsub("FMN", "fmn", demo$subject);demo$subject<-gsub("HCS", "hcs", demo$subject)
data<-right_join(demo, data, by = NULL)

data2<-read.csv(file.path(DIR, "TensorBLNAC_wm_FA_MD_BAS.csv")) #EDIT (1)
data2<-data2[,c("subject", "group", "bas_total", "basrr_total", "basfs_total", "basdrive_total", 
                "nac_mpfc_Mean_FA", "nac_mpfc_Mean_MD")]
data2<-right_join(demo, data2, by = NULL)

data3<-read.csv(file.path(DIR, "nac_mpfc_FA_MD_behavior.csv")) #EDIT (1)
data3$subject<-ifelse(grepl("FMN", data3$subject), 
                      sub("FMN", "fmn", data3$subject), sub("HCS", "hcs", data3$subject))
data2<-left_join(data2, data3, by = NULL)
data2 <- data2 %>%
  mutate(across(12:23, as.numeric))

#JHU data
data3<-read_excel(file.path(DIR, "OPAL01_JHU_primary_FA_MD_BAS.xlsx")) #EDIT (1)
data3.1<-subset(data3, select = -c(MD))
data3.1<-spread(data3.1, key = "roi", value = "FA")
colnames(data3.1)[9:15]<-paste(colnames(data3.1)[9:15], "FA", sep = "_")
data3.2<-subset(data3, select = -c(FA))
data3.2<-spread(data3.2, key = "roi", value = "MD")
colnames(data3.2)[9:15]<-paste(colnames(data3.2)[9:15], "MD", sep = "_")
data3<-full_join(data3.1, data3.2, by = NULL);rm(data3.1, data3.2)

data4<-read.csv(file.path(DIR, "OPAL01_JHU_supp.csv")) #EDIT (1)
colnames(data4)[1]<-"subject"
data4$subject<-ifelse(grepl("FMN", data4$subject), 
                      sub("FMN", "fmn", data4$subject), sub("HCS", "hcs", data4$subject))
data3<-full_join(data3, data4, by = NULL);rm(data4)

table(data$group)
# fmn hcs 
# 21  21

#######################
#Primary
#######################
MDAT2<-data2[,c(1,4,3, 9:10)]
MDAT2$group<-relevel(as.factor(MDAT2$group), ref = "HCS")
mres<-lapply(4, function(i){
  dat<-na.omit(MDAT2[,c(1:3,i)]);colnames(dat)<-c("Subject", "Group", "Age", "Outcome")
  m<-lm(Outcome ~ as.factor(Group)+scale(Age), data = dat)
  m2<-lm(scale(Outcome) ~ as.factor(Group)+scale(Age), data = dat)
  an<-car::Anova(m, type = "II")
  nd<-data.frame(Group = factor(levels(dat$Group), levels = levels(dat$Group)), Age = mean(dat$Age)); adj<-setNames(predict(m, nd), levels(dat$Group)) #EDIT (4) age-adjusted means at the sample's mean age
  
  #getting descriptive stats
  sum<-dat %>% group_by(Group) %>%
    summarise(
      N = n(),
      AdjMean = unname(adj[as.character(Group[1])]), #EDIT (4) age-adjusted mean
      AdjSD = sd(Outcome) #EDIT (4) observed SD
    )
  sum$Outcome<-colnames(MDAT2)[i]
  sum<-gather(sum, key = "Meas", value = "val", N:AdjSD)
  sum$Label<-paste(sum$Meas, sum$Group, sep = "_")
  sum<-subset(sum, select = -c(Meas, Group));sum<-spread(sum, key = "Label", value = "val")

  sum$Adj_SDP<-sqrt(((sum$N_FMN - 1)* sum$AdjSD_FMN^2 +
                       (sum$N_HCS - 1)* sum$AdjSD_HCS^2 ) /
                      (sum$N_FMN +sum$N_HCS -2))
  sum$AdjSMD<- (sum$AdjMean_FMN - sum$AdjMean_HCS) / sum$Adj_SDP
  sum<-sum[,c("Outcome", "N_FMN", "N_HCS", 
              "AdjMean_FMN","AdjMean_HCS", "AdjSD_FMN","AdjSD_HCS","AdjSMD")]
  sum$AdjMean_FMN<-round(sum$AdjMean_FMN,6);sum$AdjMean_HCS<-round(sum$AdjMean_HCS,6)
  sum$AdjSD_FMN<-round(sum$AdjSD_FMN,6);sum$AdjSD_HCS<-round(sum$AdjSD_HCS,6)
  sum$AdjSMD<-round(sum$AdjSMD,6)
  fin<-as.data.frame(
    cbind(
      sum,
      F   = round(an[1,3], 3),
      df1 = an[1,2], #group
      df2 = an[3,2], #residuals
      p   = round(an[1,4], 4),
      Beta = round(coef(m)[2],6),
      Lower = confint(m)[2,1],
      Upper = confint(m)[2,2],
      BetaS = round(coef(m2)[2],3),
      LowerS = confint(m2)[2,1],
      UpperS = confint(m2)[2,2]
    )
  )
  fin
  
})
mres2<-lapply(5, function(i){
  dat<-na.omit(MDAT2[,c(1:3,i)]);colnames(dat)<-c("Subject", "Group", "Age", "Outcome")
  m<-lm(I(Outcome*1e5) ~ as.factor(Group)+scale(Age), data = dat)
  m2<-lm(scale(Outcome) ~ as.factor(Group)+scale(Age), data = dat)
  an<-car::Anova(m, type = "II")
  nd<-data.frame(Group = factor(levels(dat$Group), levels = levels(dat$Group)), Age = mean(dat$Age)); adj<-setNames(predict(m, nd)/1e5, levels(dat$Group)) #EDIT (4) age-adjusted means at the sample's mean age, back in MD units
  
  #getting descriptive stats
  sum<-dat %>% group_by(Group) %>%
    summarise(
      N = n(),
      AdjMean = unname(adj[as.character(Group[1])]), #EDIT (4) age-adjusted mean
      AdjSD = sd(Outcome) #EDIT (4) observed SD
    )
  sum$Outcome<-colnames(MDAT2)[i]
  sum<-gather(sum, key = "Meas", value = "val", N:AdjSD)
  sum$Label<-paste(sum$Meas, sum$Group, sep = "_")
  sum<-subset(sum, select = -c(Meas, Group));sum<-spread(sum, key = "Label", value = "val")

  sum$Adj_SDP<-sqrt(((sum$N_FMN - 1)* sum$AdjSD_FMN^2 +
                       (sum$N_HCS - 1)* sum$AdjSD_HCS^2 ) /
                      (sum$N_FMN +sum$N_HCS -2))
  sum$AdjSMD<- (sum$AdjMean_FMN - sum$AdjMean_HCS) / sum$Adj_SDP
  sum<-sum[,c("Outcome", "N_FMN", "N_HCS", 
              "AdjMean_FMN","AdjMean_HCS", "AdjSD_FMN","AdjSD_HCS","AdjSMD")]

  sum$AdjMean_FMN<-round(sum$AdjMean_FMN,6);sum$AdjMean_HCS<-round(sum$AdjMean_HCS,6)
  sum$AdjSD_FMN<-round(sum$AdjSD_FMN,6);sum$AdjSD_HCS<-round(sum$AdjSD_HCS,6)
  sum$AdjSMD<-round(sum$AdjSMD,6)
  fin<-as.data.frame(
    cbind(
      sum,
      F   = round(an[1,3], 3),
      df1 = an[1,2], #group
      df2 = an[3,2], #residuals
      p   = round(an[1,4], 4),
      Beta = round(coef(m)[2],6),
      Lower = confint(m)[2,1],
      Upper = confint(m)[2,2],
      BetaS = round(coef(m2)[2],3),
      LowerS = confint(m2)[2,1],
      UpperS = confint(m2)[2,2]
    )
  )
  fin
  fin$Beta<-I(fin$Beta/1e5)
  fin$Lower<-I(fin$Lower/1e5)
  fin$Upper<-I(fin$Upper/1e5)
  fin
  
})
mres<-rbind(mres[[1]], mres2[[1]]);rm(mres2)
#all have at least 17 so keep call for correction
#correct for multiple comparisons (Holm)
mresO<-mres
mres$p.adj<-p.adjust(mres$p, method = "holm")

mres1<-mres[grepl("MD", mres$Outcome),]
mres1$AdjMean_FMN<-round(mres1$AdjMean_FMN,7)
mres1$AdjMean_HCS<-round(mres1$AdjMean_HCS,7)
mres1$AdjSD_FMN<-round(mres1$AdjSD_FMN,7)
mres1$AdjSD_HCS<-round(mres1$AdjSD_HCS,7)
mres1$Beta<-round(mres1$Beta,7)
mres1$Lower<-round(mres1$Lower,7)
mres1$Upper<-round(mres1$Upper,7)
mres1$Mean_FMN<-paste0(mres1$AdjMean_FMN*10000,"e-4", " (", mres1$AdjSD_FMN, ")")
mres1$Mean_HCS<-paste0(mres1$AdjMean_HCS*10000,"e-4",  " (", mres1$AdjSD_HCS,")")
mres1$Est<-paste0(mres1$Beta," (", mres1$Lower,", ", mres1$Upper,")")
mres1$Est_Std<-paste0(mres1$BetaS, " (", round(mres1$LowerS,3),", ", round(mres1$UpperS,3), ")")

mres2<-mres[grepl("FA", mres$Outcome),]
mres2$AdjMean_FMN<-round(mres2$AdjMean_FMN,3);mres2$AdjMean_HCS<-round(mres2$AdjMean_HCS,3)
mres2$AdjSD_FMN<-round(mres2$AdjSD_FMN,3);mres2$AdjSD_HCS<-round(mres2$AdjSD_HCS,3)
mres2$Beta<-round(mres2$Beta,3)
mres2$Mean_FMN<-paste0(mres2$AdjMean_FMN, " (", mres2$AdjSD_FMN, ")")
mres2$Mean_HCS<-paste0(mres2$AdjMean_HCS, " (", mres2$AdjSD_HCS, ")")
mres2$Est<-paste0(mres2$Beta, " (", round(mres2$Lower,3),", ", round(mres2$Upper,3), ")")
mres2$Est_Std<-paste0(round(mres2$BetaS,3), " (", round(mres2$LowerS,3),", ", round(mres2$UpperS,3), ")")

mres<-rbind(mres1, mres2);rm(mres1, mres2)
mres$AdjSMD<-round(mres$AdjSMD,3)
mresC<-mres[,c("Outcome", "N_FMN", "N_HCS","AdjMean_FMN", "AdjMean_HCS", "Mean_FMN", "Mean_HCS", "AdjSMD","Est","Est_Std", "F", "df1", "df2", "p", "p.adj")]

Prim1<-list(mresO, mres, mresC);rm(mresO, mres, mresC)

#######################
#Exploratory
#######################
MDAT<-data[,c(1,4,3, 9:22)]
MDAT$group<-relevel(as.factor(MDAT$group), ref = "HCS")
mres<-lapply(4:10, function(i){
  dat<-na.omit(MDAT[,c(1:3,i)]);colnames(dat)<-c("Subject", "Group", "Age", "Outcome")
  m<-lm(Outcome ~ as.factor(Group)+scale(Age), data = dat)
  m2<-lm(scale(Outcome) ~ as.factor(Group)+scale(Age), data = dat)
  an<-car::Anova(m, type = "II")
  nd<-data.frame(Group = factor(levels(dat$Group), levels = levels(dat$Group)), Age = mean(dat$Age)); adj<-setNames(predict(m, nd), levels(dat$Group)) #EDIT (4) age-adjusted means at the sample's mean age
  
  #getting descriptive stats
  sum<-dat %>% group_by(Group) %>%
    summarise(
      N = n(),
      AdjMean = unname(adj[as.character(Group[1])]), #EDIT (4) age-adjusted mean
      AdjSD = sd(Outcome) #EDIT (4) observed SD
    )
  sum$Outcome<-colnames(MDAT)[i]
  sum<-gather(sum, key = "Meas", value = "val", N:AdjSD)
  sum$Label<-paste(sum$Meas, sum$Group, sep = "_")
  sum<-subset(sum, select = -c(Meas, Group));sum<-spread(sum, key = "Label", value = "val")

  sum$Adj_SDP<-sqrt(((sum$N_FMN - 1)* sum$AdjSD_FMN^2 +
                       (sum$N_HCS - 1)* sum$AdjSD_HCS^2 ) /
                      (sum$N_FMN +sum$N_HCS -2))
  sum$AdjSMD<- (sum$AdjMean_FMN - sum$AdjMean_HCS) / sum$Adj_SDP
  sum<-sum[,c("Outcome", "N_FMN", "N_HCS", 
              "AdjMean_FMN","AdjMean_HCS", "AdjSD_FMN","AdjSD_HCS","AdjSMD")]
  sum$AdjMean_FMN<-round(sum$AdjMean_FMN,6);sum$AdjMean_HCS<-round(sum$AdjMean_HCS,6)
  sum$AdjSD_FMN<-round(sum$AdjSD_FMN,6);sum$AdjSD_HCS<-round(sum$AdjSD_HCS,6)
  sum$AdjSMD<-round(sum$AdjSMD,6)
  fin<-as.data.frame(
    cbind(
      sum,
      F   = round(an[1,3], 3),
      df1 = an[1,2], #group
      df2 = an[3,2], #residuals
      p   = round(an[1,4], 4),
      Beta = round(coef(m)[2],6),
      Lower = confint(m)[2,1],
      Upper = confint(m)[2,2],
      BetaS = round(coef(m2)[2],3),
      LowerS = confint(m2)[2,1],
      UpperS = confint(m2)[2,2]
    )
  )
  fin
  
})
mres<-do.call("rbind", mres)

mres2<-lapply(11:17, function(i){
  dat<-na.omit(MDAT[,c(1:3,i)]);colnames(dat)<-c("Subject", "Group", "Age", "Outcome")
  m<-lm(I(Outcome*1e5) ~ as.factor(Group)+scale(Age), data = dat)
  m2<-lm(scale(Outcome) ~ as.factor(Group)+scale(Age), data = dat)
  an<-car::Anova(m, type = "II")
  nd<-data.frame(Group = factor(levels(dat$Group), levels = levels(dat$Group)), Age = mean(dat$Age)); adj<-setNames(predict(m, nd)/1e5, levels(dat$Group)) #EDIT (4) age-adjusted means at the sample's mean age, back in MD units
  
  #getting descriptive stats
  sum<-dat %>% group_by(Group) %>%
    summarise(
      N = n(),
      # RawMean = mean(Outcome),
      # RawSD = sd(Outcome),
      AdjMean = unname(adj[as.character(Group[1])]), #EDIT (4) age-adjusted mean
      AdjSD = sd(Outcome) #EDIT (4) observed SD
    )
  sum$Outcome<-colnames(MDAT)[i]
  sum<-gather(sum, key = "Meas", value = "val", N:AdjSD)
  sum$Label<-paste(sum$Meas, sum$Group, sep = "_")
  sum<-subset(sum, select = -c(Meas, Group));sum<-spread(sum, key = "Label", value = "val")

  sum$Adj_SDP<-sqrt(((sum$N_FMN - 1)* sum$AdjSD_FMN^2 +
                       (sum$N_HCS - 1)* sum$AdjSD_HCS^2 ) /
                      (sum$N_FMN +sum$N_HCS -2))
  sum$AdjSMD<- (sum$AdjMean_FMN - sum$AdjMean_HCS) / sum$Adj_SDP
  sum<-sum[,c("Outcome", "N_FMN", "N_HCS", 
              "AdjMean_FMN","AdjMean_HCS", "AdjSD_FMN","AdjSD_HCS","AdjSMD")]

  sum$AdjMean_FMN<-round(sum$AdjMean_FMN,6);sum$AdjMean_HCS<-round(sum$AdjMean_HCS,6)
  sum$AdjSD_FMN<-round(sum$AdjSD_FMN,6);sum$AdjSD_HCS<-round(sum$AdjSD_HCS,6)
  sum$AdjSMD<-round(sum$AdjSMD,6)
  fin<-as.data.frame(
    cbind(
      sum,
      F   = round(an[1,3], 3),
      df1 = an[1,2], #group
      df2 = an[3,2], #residuals
      p   = round(an[1,4], 4),
      Beta = round(coef(m)[2],6),
      Lower = confint(m)[2,1],
      Upper = confint(m)[2,2],
      BetaS = round(coef(m2)[2],3),
      LowerS = confint(m2)[2,1],
      UpperS = confint(m2)[2,2]
    )
  )
  fin
  fin$Beta<-I(fin$Beta/1e5)
  fin$Lower<-I(fin$Lower/1e5)
  fin$Upper<-I(fin$Upper/1e5)
  fin
  
})
mres2<-do.call("rbind", mres2)
mres<-rbind(mres, mres2);rm(mres2)
#all have at least 17 so keep call for correction
mresO<-mres
#corrected within FA and within MD
mres1<-mres[grepl("MD", mres$Outcome),];mres1$p.adj<-p.adjust(mres1$p, method = "BH")
mres1$AdjMean_FMN<-round(mres1$AdjMean_FMN,7)
mres1$AdjMean_HCS<-round(mres1$AdjMean_HCS,7)
mres1$AdjSD_FMN<-round(mres1$AdjSD_FMN,7)
mres1$AdjSD_HCS<-round(mres1$AdjSD_HCS,7)
mres1$Beta<-round(mres1$Beta,7)
mres1$Lower<-round(mres1$Lower,7)
mres1$Upper<-round(mres1$Upper,7)
mres1$Mean_FMN<-paste0(mres1$AdjMean_FMN*10000,"e-4", " (", mres1$AdjSD_FMN, ")")
mres1$Mean_HCS<-paste0(mres1$AdjMean_HCS*10000,"e-4",  " (", mres1$AdjSD_HCS,")")
mres1$Est<-paste0(mres1$Beta," (", mres1$Lower,", ", mres1$Upper,")")
mres1$Est_Std<-paste0(mres1$BetaS, " (", round(mres1$LowerS,3),", ", round(mres1$UpperS,3), ")")

mres2<-mres[grepl("FA", mres$Outcome),];mres2$p.adj<-p.adjust(mres2$p, method = "BH")
mres2$AdjMean_FMN<-round(mres2$AdjMean_FMN,3);mres2$AdjMean_HCS<-round(mres2$AdjMean_HCS,3)
mres2$AdjSD_FMN<-round(mres2$AdjSD_FMN,3);mres2$AdjSD_HCS<-round(mres2$AdjSD_HCS,3)
mres2$Beta<-round(mres2$Beta,3)
mres2$Mean_FMN<-paste0(mres2$AdjMean_FMN, " (", mres2$AdjSD_FMN, ")")
mres2$Mean_HCS<-paste0(mres2$AdjMean_HCS, " (", mres2$AdjSD_HCS, ")")
mres2$Est<-paste0(mres2$Beta, " (", round(mres2$Lower,3),", ", round(mres2$Upper,3), ")")
mres2$Est_Std<-paste0(round(mres2$BetaS,3), " (", round(mres2$LowerS,3),", ", round(mres2$UpperS,3), ")")


mres<-rbind(mres1, mres2);rm(mres1, mres2)
mres$AdjSMD<-round(mres$AdjSMD,3)
mresC<-mres[,c("Outcome", "N_FMN", "N_HCS","AdjMean_FMN", "AdjMean_HCS", "Mean_FMN", "Mean_HCS", "AdjSMD","Est","Est_Std", "F", "df1", "df2", "p", "p.adj")]

Supp<-list(mresO, mres, mresC);rm(mresO, mres, mresC)

#######################
#JHU
#######################
data3$group<-relevel(as.factor(data3$group), ref = "HCS")
mres<-lapply(c(9:15,27,29,31,33,35,37), function(i){
  #print(i)
  dat<-na.omit(data3[,c(1:3,i)]);colnames(dat)<-c("Subject", "Group", "Age", "Outcome")
  m<-lm(Outcome ~ as.factor(Group)+scale(Age), data = dat)
  m2<-lm(scale(Outcome) ~ as.factor(Group)+scale(Age), data = dat)
  an<-car::Anova(m, type = "II")
  nd<-data.frame(Group = factor(levels(dat$Group), levels = levels(dat$Group)), Age = mean(dat$Age)); adj<-setNames(predict(m, nd), levels(dat$Group)) #EDIT (4) age-adjusted means at the sample's mean age
  
  #getting descriptive stats
  sum<-dat %>% group_by(Group) %>%
    summarise(
      N = n(),
      AdjMean = unname(adj[as.character(Group[1])]), #EDIT (4) age-adjusted mean
      AdjSD = sd(Outcome) #EDIT (4) observed SD
    )
  sum$Outcome<-colnames(data3)[i]
  sum<-gather(sum, key = "Meas", value = "val", N:AdjSD)
  sum$Label<-paste(sum$Meas, sum$Group, sep = "_")
  sum<-subset(sum, select = -c(Meas, Group));sum<-spread(sum, key = "Label", value = "val")

  sum$Adj_SDP<-sqrt(((sum$N_FMN - 1)* sum$AdjSD_FMN^2 +
                       (sum$N_HCS - 1)* sum$AdjSD_HCS^2 ) /
                      (sum$N_FMN +sum$N_HCS -2))
  sum$AdjSMD<- (sum$AdjMean_FMN - sum$AdjMean_HCS) / sum$Adj_SDP
  sum<-sum[,c("Outcome", "N_FMN", "N_HCS", 
              "AdjMean_FMN","AdjMean_HCS", "AdjSD_FMN","AdjSD_HCS","AdjSMD")]
  sum$AdjMean_FMN<-round(sum$AdjMean_FMN,6);sum$AdjMean_HCS<-round(sum$AdjMean_HCS,6)
  sum$AdjSD_FMN<-round(sum$AdjSD_FMN,6);sum$AdjSD_HCS<-round(sum$AdjSD_HCS,6)
  sum$AdjSMD<-round(sum$AdjSMD,6)
  fin<-as.data.frame(
    cbind(
      sum,
      F   = round(an[1,3], 3),
      df1 = an[1,2], #group
      df2 = an[3,2], #residuals
      p   = round(an[1,4], 4),
      Beta = round(coef(m)[2],6),
      Lower = confint(m)[2,1],
      Upper = confint(m)[2,2],
      BetaS = round(coef(m2)[2],3),
      LowerS = confint(m2)[2,1],
      UpperS = confint(m2)[2,2]
    )
  )
  fin
  
})
mres<-do.call("rbind", mres)

mres2<-lapply(c(16:22,28,30,32,34,36,38), function(i){
  dat<-na.omit(data3[,c(1:3,i)]);colnames(dat)<-c("Subject", "Group", "Age", "Outcome")
  m<-lm(I(Outcome*1e5) ~ as.factor(Group)+scale(Age), data = dat)
  m2<-lm(scale(Outcome) ~ as.factor(Group)+scale(Age), data = dat)
  an<-car::Anova(m, type = "II")
  nd<-data.frame(Group = factor(levels(dat$Group), levels = levels(dat$Group)), Age = mean(dat$Age)); adj<-setNames(predict(m, nd)/1e5, levels(dat$Group)) #EDIT (4) age-adjusted means at the sample's mean age, back in MD units
  
  #getting descriptive stats
  sum<-dat %>% group_by(Group) %>%
    summarise(
      N = n(),
      AdjMean = unname(adj[as.character(Group[1])]), #EDIT (4) age-adjusted mean
      AdjSD = sd(Outcome) #EDIT (4) observed SD
    )
  sum$Outcome<-colnames(data3)[i]
  sum<-gather(sum, key = "Meas", value = "val", N:AdjSD)
  sum$Label<-paste(sum$Meas, sum$Group, sep = "_")
  sum<-subset(sum, select = -c(Meas, Group));sum<-spread(sum, key = "Label", value = "val")
  #getting SMD
  sum$Adj_SDP<-sqrt(((sum$N_FMN - 1)* sum$AdjSD_FMN^2 +
                       (sum$N_HCS - 1)* sum$AdjSD_HCS^2 ) /
                      (sum$N_FMN +sum$N_HCS -2))
  sum$AdjSMD<- (sum$AdjMean_FMN - sum$AdjMean_HCS) / sum$Adj_SDP
  sum<-sum[,c("Outcome", "N_FMN", "N_HCS", 
              "AdjMean_FMN","AdjMean_HCS", "AdjSD_FMN","AdjSD_HCS","AdjSMD")]
  sum$AdjMean_FMN<-round(sum$AdjMean_FMN,6);sum$AdjMean_HCS<-round(sum$AdjMean_HCS,6)
  sum$AdjSD_FMN<-round(sum$AdjSD_FMN,6);sum$AdjSD_HCS<-round(sum$AdjSD_HCS,6)
  sum$AdjSMD<-round(sum$AdjSMD,6)
  fin<-as.data.frame(
    cbind(
      sum,
      F   = round(an[1,3], 3),
      df1 = an[1,2], #group
      df2 = an[3,2], #residuals
      p   = round(an[1,4], 4),
      Beta = round(coef(m)[2],6),
      Lower = confint(m)[2,1],
      Upper = confint(m)[2,2],
      BetaS = round(coef(m2)[2],3),
      LowerS = confint(m2)[2,1],
      UpperS = confint(m2)[2,2]
    )
  )
  fin
  fin$Beta<-I(fin$Beta/1e5)
  fin$Lower<-I(fin$Lower/1e5)
  fin$Upper<-I(fin$Upper/1e5)
  fin
  
})
mres2<-do.call("rbind", mres2)
mres<-rbind(mres, mres2);rm(mres2)
#corrected within FA and within MD
mres1<-mres[grepl("MD", mres$Outcome),]
mres1$AdjMean_FMN<-round(mres1$AdjMean_FMN,7)
mres1$AdjMean_HCS<-round(mres1$AdjMean_HCS,7)
mres1$AdjSD_FMN<-round(mres1$AdjSD_FMN,7)
mres1$AdjSD_HCS<-round(mres1$AdjSD_HCS,7)
mres1$Beta<-round(mres1$Beta,7)
mres1$Lower<-round(mres1$Lower,7)
mres1$Upper<-round(mres1$Upper,7)
mres1$Mean_FMN<-paste0(mres1$AdjMean_FMN*10000,"e-4", " (", mres1$AdjSD_FMN, ")")
mres1$Mean_HCS<-paste0(mres1$AdjMean_HCS*10000,"e-4",  " (", mres1$AdjSD_HCS,")")
mres1$Est<-paste0(mres1$Beta," (", mres1$Lower,", ", mres1$Upper,")")
mres1$Est_Std<-paste0(mres1$BetaS, " (", round(mres1$LowerS,3),", ", round(mres1$UpperS,3), ")")

mres2<-mres[grepl("FA", mres$Outcome),]
mres2$AdjMean_FMN<-round(mres2$AdjMean_FMN,3);mres2$AdjMean_HCS<-round(mres2$AdjMean_HCS,3)
mres2$AdjSD_FMN<-round(mres2$AdjSD_FMN,3);mres2$AdjSD_HCS<-round(mres2$AdjSD_HCS,3)
mres2$Beta<-round(mres2$Beta,3)
mres2$Mean_FMN<-paste0(mres2$AdjMean_FMN, " (", mres2$AdjSD_FMN, ")")
mres2$Mean_HCS<-paste0(mres2$AdjMean_HCS, " (", mres2$AdjSD_HCS, ")")
mres2$Est<-paste0(mres2$Beta, " (", round(mres2$Lower,3),", ", round(mres2$Upper,3), ")")
mres2$Est_Std<-paste0(round(mres2$BetaS,3), " (", round(mres2$LowerS,3),", ", round(mres2$UpperS,3), ")")
mresC.JHU<-rbind(mres1,mres2);rm(mres1, mres2)


##############################################
#H2 Primary (BAS_RR correlations)
##############################################
Sub<-data2[data$group == "FMN",]
#Partial Correlations for bas_rr
cres<-lapply(6, function(i){
  res<-lapply(9:10, function(j){
    dat<-na.omit(Sub[,c(1,3,i,j)])
    colnames(dat)<-c("Subject", "Age", "Bas", "Variable")
    c<-pcor.test(dat$Bas, dat$Variable, dat$Age)
    cres<-cbind(Variable1 = colnames(Sub)[i],
                Variable2 = colnames(Sub)[j],
                as.data.frame(c))
  })
  comb<-do.call("rbind", res)
  comb
})
cres<-do.call("rbind", cres)
cres$adj.p.value<-p.adjust(cres$p.value, method = "holm")

#running rest of bas
cres2<-lapply(c(5,7:8), function(i){
  res<-lapply(9:10, function(j){
    dat<-na.omit(Sub[,c(1,3,i,j)])
    colnames(dat)<-c("Subject", "Age", "Bas", "Variable")
    c<-pcor.test(dat$Bas, dat$Variable, dat$Age)
    cres<-cbind(Variable1 = colnames(Sub)[i],
                Variable2 = colnames(Sub)[j],
                as.data.frame(c))
  })
  comb<-do.call("rbind", res)
  comb
})
cres2<-do.call("rbind", cres2)
cres2$adj.p.value<-p.adjust(cres2$p.value, method = "BH")


#running shapiro-wilk tests (checking for normality)
swtest<-lapply(5:10, function(j){
  dat<-Sub[,c(j)]; colnames(dat)<-c("Var")
  x<-shapiro.test(dat$Var)
  res<-as.data.frame(cbind(Variable = colnames(Sub)[j],
                           Pvalue = x[["p.value"]]))
})
swtest<-do.call("rbind", swtest)
swtest$Pvalue<-round(as.numeric(swtest$Pvalue),3)
swtestpassed<-swtest[swtest$Pvalue >0.05,] #all passed, can just use pearson

#behavior correlations 
cres3<-lapply(9, function(i){
  res<-lapply(12:23, function(j){
    dat<-na.omit(Sub[,c(1,3,i,j)])
    colnames(dat)<-c("Subject", "Age", "Var1", "Var2")
    c<-cor.test(dat$Var1, dat$Var2)
    cres<-as.data.frame(cbind(Variable1 = colnames(Sub)[i],
                Variable2 = colnames(Sub)[j],
                CorEst = c[["estimate"]][["cor"]],
                Pvalue = c[["p.value"]]))
  })
  comb<-do.call("rbind", res)
  comb
})
cres3<-do.call("rbind", cres3)
cres3$CorEst<-round(as.numeric(cres3$CorEst),3)
cres3$Pvalue<-as.numeric(cres3$Pvalue)

cres3.2<-lapply(9, function(i){
  res<-lapply(12:23, function(j){
    dat<-na.omit(Sub[,c(1,3,i,j)])
    colnames(dat)<-c("Subject", "Age", "Var1", "Var2")
    c<-cor.test(dat$Var1, dat$Var2, method = "spearman", exact = F)
    cres<-as.data.frame(cbind(Variable1 = colnames(Sub)[i],
                              Variable2 = colnames(Sub)[j],
                              CorEst = c[["estimate"]],
                              Pvalue = c[["p.value"]]))
  })
  comb<-do.call("rbind", res)
  comb
})
cres3.2<-do.call("rbind", cres3.2)
cres3.2$CorEst<-round(as.numeric(cres3.2$CorEst),3)
cres3.2$Pvalue<-as.numeric(cres3.2$Pvalue)

#running shapiro-wilk tests (checking for normality)
swtest<-lapply(12:23, function(j){
  dat<-Sub[,c(j)]; colnames(dat)<-c("Var")
  x<-shapiro.test(dat$Var)
  res<-as.data.frame(cbind(Variable = colnames(Sub)[j],
                           Pvalue = x[["p.value"]]))
})
swtest<-do.call("rbind", swtest)
swtest$Pvalue<-round(as.numeric(swtest$Pvalue),3)
swtestpassed<-swtest[swtest$Pvalue >0.05,] #all passed, can just use pearson
cres3<-cres3[cres3$Variable2 %in% swtestpassed$Variable,]
cres3.2<-cres3.2[!cres3.2$Variable2 %in% swtestpassed$Variable,]
cres3<-rbind(cres3, cres3.2);rm(cres3.2)
cres3$Test<-ifelse(cres3$Variable2 %in% swtestpassed$Variable, "Pearson", "Spearman")

write.xlsx(cres3, file.path(OUT, "BehavCorr.xlsx")) #EDIT (3)


#######################
#Saving Results
#######################

datalines<-list(
  "H1 Res" = Prim1[[3]],
  "H2" = cres,
  "H1 Exploratory" = Supp[[3]],
  "H2 Exploratory" = cres2,
  "JHU" = mresC.JHU
)

write.xlsx(datalines, file.path(OUT, "ModelRes.xlsx")) #EDIT (3)
for (nm in names(datalines)) { cat("\n==========", nm, "==========\n"); print(datalines[[nm]], row.names = FALSE) } #EDIT (6) print the model tables (same as ModelRes.xlsx)
cat("\n========== Behavior correlations (BehavCorr.xlsx) ==========\n"); print(cres3, row.names = FALSE) #EDIT (6)



PDAT<-data2[,c("subject","group","age", "nac_mpfc_Mean_FA", "nac_mpfc_Mean_MD")]
table(PDAT$group)
PDAT$group<-relevel(as.factor(PDAT$group), ref = "HCS")
m<-lm(nac_mpfc_Mean_FA ~ as.factor(group)+scale(age), data = PDAT)
PDAT$FA_ageadj<-PDAT$nac_mpfc_Mean_FA - coef(m)["scale(age)"]*as.numeric(scale(PDAT$age)) #EDIT (5) observed FA moved to the sample's mean age (was fitted values)

FAp<-ggplot(PDAT, aes(x = group, y = nac_mpfc_Mean_FA, group = group, color = group))+
  geom_hline(yintercept = c(0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52), colour = "grey92", linewidth = 0.5)+ #EDIT (7) gridlines only at the labeled y values
  geom_boxplot(width = 0.4)+
  geom_jitter(width = 0.2,shape = 17, size = 2, alpha = 0.5)+
  labs(y = "Nac-MPFC FA", x = "", title = " ")+
  scale_y_continuous(breaks = c(0.28,0.30,0.32, 0.34, 0.36,0.38,0.40,0.42, 0.44, 0.46, 0.48, 0.50, 0.52),
                     labels = c("0.28", "","0.32", "", "0.36","","0.40","", "0.44", "", "0.48", "", "0.52"),
                     limits = c(0.28, 0.52))+
  scale_color_manual(values = c("darkorange2", "darkorchid4"))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    title = element_text(size = 20, face = "bold"),
    panel.grid.major.y = element_blank(), #EDIT (7)
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  );FAp

FAp2<-ggplot(PDAT, aes(x = group, y = FA_ageadj, group = group, color = group))+
  geom_hline(yintercept = c(0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52), colour = "grey92", linewidth = 0.5)+ #EDIT (7) gridlines only at the labeled y values
  geom_boxplot(width = 0.4)+
  geom_jitter(width = 0.2,shape = 17,  size = 2, alpha = 0.5)+
  labs(y = "Age Adjusted \nNac-MPFC FA", x = "", title = " ")+
  scale_y_continuous(breaks = c(0.28,0.30,0.32, 0.34, 0.36,0.38,0.40,0.42, 0.44, 0.46, 0.48, 0.50, 0.52), #EDIT (5) same axis as the raw FA graph
                     labels = c("0.28", "","0.32", "", "0.36","","0.40","", "0.44", "", "0.48", "", "0.52"), #EDIT (5)
                     limits = c(0.28, 0.52))+ #EDIT (5)
  scale_color_manual(values = c("darkorange2", "darkorchid4"))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    title = element_text(size = 20, face = "bold"),
    panel.grid.major.y = element_blank(), #EDIT (7)
    panel.grid.minor.y = element_blank(), #EDIT (7)
    legend.position = "none"
  );FAp2


PDAT$MDrev<-PDAT$nac_mpfc_Mean_MD*10000

#7.38 to 8.84
MDp<-ggplot(PDAT, aes(x = group, y = MDrev, group = group, color = group))+
  geom_hline(yintercept = c(7.2, 7.6, 8.0, 8.4, 8.8, 9.2), colour = "grey92", linewidth = 0.5)+ #EDIT (7) gridlines only at the labeled y values
  geom_boxplot(width = 0.4)+
  geom_jitter(width = 0.2, size = 2, alpha = 0.5)+
  labs(y = "Nac-MPFC FA \nMean Diffusivity", x = "", title = "C")+
  scale_y_continuous(breaks = c(7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0,9.2),
                     labels = c(expression("7.2 x10"^~"-4"),"",
                                expression("7.6 x10"^~"-4"), "",
                                expression("8.0 x10"^~"-4"),"",
                                expression("8.4 x10"^~"-4"),"",
                                expression("8.8 x10"^~"-4"), "",
                                expression("9.2 x10"^~"-4")),
                     limits = c(7.2, 9.2))+
  scale_color_manual(values = c("darkorange2", "darkorchid4"))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    title = element_text(size = 20, face = "bold"),
    panel.grid.major.y = element_blank(), #EDIT (7)
    panel.grid.minor.y = element_blank(), #EDIT (7)
    legend.position = "none"
  );MDp

comboplot<-ggarrange(FAp, MDp, nrow = 1);comboplot
png(file.path(OUT, "DTI_comboplot_FA_MD.png"),width = 900, height = 400,res = 150) #EDIT (3)
print(comboplot) #EDIT (3)
dev.off() #EDIT (3)

altcomboplot<-ggarrange(FAp, FAp2, nrow = 1);altcomboplot
png(file.path(OUT, "DTI_comboplot_FA_FAadj.png"),width = 900, height = 400,res = 140) #EDIT (3)
print(altcomboplot) #EDIT (3)
dev.off() #EDIT (3)

png(file.path(OUT, "DTI_FAonly.png"),width = 600, height = 400,res = 140) #EDIT (3)
print(FAp) #EDIT (3)
dev.off() #EDIT (3)

png(file.path(OUT, "DTI_FAAdjonly.png"),width = 600, height = 400,res = 140) #EDIT (3)
print(FAp2) #EDIT (3)
dev.off() #EDIT (3)

print(FAp); print(FAp2); print(MDp); print(comboplot); print(altcomboplot) #EDIT (6) show the plots when the script is source()d
message("\nSaved to ", normalizePath(OUT), ":\n  ", paste(list.files(OUT), collapse = "\n  ")) #EDIT (6)
