################################################################################
# SCRIPT FOR PLOTTING AND ANALYZING DPRIME FROM EXPERIMENT III
################################################################################

# CTRL + ALT + R

################################################################################
# LIBRARY LOAD
################################################################################

library(IDPmisc)
library(nlme)
library(multcomp)
library(xtable)
library(car)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(doBy)
library(pracma)
library(lattice)
library(emmeans)
library(phia)
library(lme4)
library(lmerTest)
library(leaps)
library(dplyr)
library(stringr)
library(report)

################################################################################
# FUNCTION
################################################################################

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, 
                      conf.interval=.95) {
  library(doBy)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), 
                              sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean", sep="")] = measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd", sep="")] = "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="")] = "N"
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

################################################################################
# REFERENCING FOR PLOTTING
################################################################################

root_path = "/home/link/Documents/thèse_onera/"
article_path = "/home/link/Documents/thèse_onera/articles_alex/"
article_target = "article_1_PA/PloSOne/maxi_dossier/Exp_III/performance_stats"
slides_path = "/home/link/Documents/thèse_onera/diapos_phd_thesis/"
exp_path = "/home/link/Documents/thèse_onera/experimentation/"
slides_path_acc = "images/EEG/info_content"
slides_path_final = str_c(slides_path, slides_path_acc, sep="")

multcomp_adjusted_type = "fdr"
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")

################################################################################
# MODELE d' ~ MPPO
################################################################################

setwd(exp_path)

dprime_Exp_III_mppo = read.table(
  "psychoac_mppo_Exp_III.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_III_mppo = dprime_Exp_III_mppo[which(dprime_Exp_III_mppo$Sujet!=8),]

cond_mppo = c("16fpo","32fpo","64fpo")
dprime_Exp_III_mppo$MPPO = factor(dprime_Exp_III_mppo$MPPO,levels=cond_mppo)

qqnorm(dprime_Exp_III_mppo[,"dprime"], main = "QQNORM FOR Dprime MPPO")
qqline(dprime_Exp_III_mppo[,"dprime"], main = "QQNORM FOR Dprime MPPO")

mixmodel = lme(dprime~MPPO, data=dprime_Exp_III_mppo, 
               random=~1|Sujet, method="ML")
summary(mixmodel)
anova(mixmodel)

mcp_mppo = glht(mixmodel, linfct=mcp(MPPO="Tukey"))
summary(mcp_mppo,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_mppo)$test
mtests_mppo = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_mppo$alternativ, 
  less = paste("Pr(<", ifelse(mcp_mppo$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_mppo$df == 0, "z", "t"), ")", sep = ""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo$df == 0, "z", "t"), "|)", sep =""))
colnames(mtests_mppo) = c("Estimate", "Std. Error", 
                          ifelse(mcp_mppo$df ==0, "z value", "t value"), pname)

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_mppo)

################################################################################
# MODELE d' ~ MITI
################################################################################

setwd(exp_path)

dprime_Exp_III_miti = read.table(
  "psychoac_miti_Exp_III.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_III_miti = dprime_Exp_III_miti[which(dprime_Exp_III_miti$Sujet!=8),]

cond_miti = c("2300ms","1100ms","300ms")
dprime_Exp_III_miti$MITI = factor(dprime_Exp_III_miti$MITI, levels=cond_miti)

qqnorm(dprime_Exp_III_miti[,"dprime"], main = "QQNORM FOR Dprime MITI")
qqline(dprime_Exp_III_miti[,"dprime"], main = "QQNORM FOR Dprime MITI")

mixmodel = lme(dprime~MITI, data=dprime_Exp_III_miti, 
               random=~1|Sujet, method="ML")
summary(mixmodel)
anova(mixmodel)

mcp_miti = glht(mixmodel, linfct=mcp(MITI="Tukey"))
summary(mcp_miti,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_miti)$test
mtests_miti = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_miti$alternativ,
  less = paste("Pr(<", ifelse(mcp_miti$df ==0, "z", "t"), ")", sep = ""), 
  greater = paste("Pr(>", ifelse(mcp_miti$df == 0, "z", "t"), ")", sep = ""), 
  two.sided = paste("Pr(>|", ifelse(mcp_miti$df == 0, "z", "t"), "|)", sep = ""))                                                                   
colnames(mtests_miti) = c("Estimate", "Std. Error", 
                         ifelse(mcp_miti$df ==0, "z value", "t value"), pname)
report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_miti)

################################################################################
# MODELE d' ~ TARGET RATE
################################################################################

setwd(exp_path)

dprime_Exp_III_tra = read.table(
  "psychoac_targetrate_Exp_III.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_III_tra = dprime_Exp_III_tra[which(dprime_Exp_III_tra$Sujet!=8),]

cond_tra = c("1Hz","2Hz","5Hz")
dprime_Exp_III_tra$TargetRate = factor(dprime_Exp_III_tra$TargetRate, 
                                      levels=cond_tra)

qqnorm(dprime_Exp_III_tra[,"dprime"], main = "QQNORM FOR Dprime TargetRate")
qqline(dprime_Exp_III_tra[,"dprime"], main = "QQNORM FOR Dprime TargetRate")

mixmodel = lme(dprime~TargetRate, data=dprime_Exp_III_tra, 
               random=~1|Sujet, method="ML")
summary(mixmodel)
anova(mixmodel)

mcp_tra = glht(mixmodel, linfct=mcp(TargetRate="Tukey"))
summary(mcp_tra,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_tra)$test
mtests_tra = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""), 
  greater = paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep = ""), 
  two.sided = paste("Pr(>|", ifelse(mcp_tra$df == 0, "z", "t"), "|)", sep = ""))
colnames(mtests_tra) = c("Estimate", "Std. Error", 
                         ifelse(mcp_tra$df ==0, "z value", "t value"), pname)
report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_tra)

################################################################################
# MODELE d' ~ UNCERTAINTY
################################################################################

setwd(exp_path)

dprime_Exp_III_uncertainty = read.table(
  "psychoac_uncertainty_Exp_III.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_III_uncertainty = dprime_Exp_III_uncertainty[which(
  dprime_Exp_III_uncertainty$Sujet!=8),]

cond_uncertainty = c("84","110","123","169","221","246","339","442","492")
dprime_Exp_III_uncertainty$Uncertainty = factor(
  dprime_Exp_III_uncertainty$Uncertainty,levels=cond_uncertainty)

qqnorm(dprime_Exp_III_uncertainty[,"dprime"], 
       main = "QQNORM FOR Dprime Uncertainty")
qqline(dprime_Exp_III_uncertainty[,"dprime"], 
       main = "QQNORM FOR Dprime Uncertainty")

mixmodel = lme(dprime~Uncertainty, data=dprime_Exp_III_uncertainty, 
               random=~1|Sujet, method="ML")
summary(mixmodel)
anova(mixmodel)

mcp_uncertainty = glht(mixmodel, linfct=mcp(Uncertainty="Tukey"))
summary(mcp_uncertainty,test=adjusted(type="bonferroni"))
pq = summary(mcp_uncertainty)$test
mtests_uncertainty = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_uncertainty$alternativ, 
  less = paste("Pr(<", ifelse(
    mcp_uncertainty$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(
    mcp_uncertainty$df == 0, "z", "t"), ")", sep = ""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_uncertainty$df == 0, "z", "t"), "|)", sep =""))
colnames(mtests_uncertainty) = c("Estimate", "Std. Error", ifelse(
  mcp_uncertainty$df ==0, "z value", "t value"), pname)

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_uncertainty)

################################################################################
# MODELE d' ~ MPPO * MITI
################################################################################

setwd(exp_path)

dprime_Exp_III_mppo_miti = read.table(
  "psychoac_mppo_miti_Exp_III.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_III_mppo_miti = dprime_Exp_III_mppo_miti[which(
  dprime_Exp_III_mppo_miti$Sujet!=8),]

cond_mppo = c("16fpo","32fpo","64fpo")
dprime_Exp_III_mppo_miti$MPPO = factor(
  dprime_Exp_III_mppo_miti$MPPO, levels=cond_mppo)

cond_miti = c("2300ms","1100ms","300ms")
dprime_Exp_III_mppo_miti$MITI = factor(
  dprime_Exp_III_mppo_miti$MITI, levels=cond_miti)

qqnorm(dprime_Exp_III_mppo_miti[,"dprime"], main = "QQNORM FOR Dprime")
qqline(dprime_Exp_III_mppo_miti[,"dprime"], main = "QQNORM FOR Dprime")

mixmodel = lme(dprime ~ MPPO * MITI, data=dprime_Exp_III_mppo_miti,
               random=~1|Sujet, method= "ML")

summary(mixmodel)
anova(mixmodel)

dprime_Exp_III_mppo_miti$res = interaction(dprime_Exp_III_mppo_miti$MPPO, 
                                        dprime_Exp_III_mppo_miti$MITI)
dprime_Exp_III_mppo_miti$res = factor(dprime_Exp_III_mppo_miti$res)

res_mixmodel = lme(dprime ~ res, data=dprime_Exp_III_mppo_miti, random=~1|Sujet, 
                   method= "ML")

mcp_mppo_miti = glht(res_mixmodel, linfct=mcp(res="Tukey"))
summary(mcp_mppo_miti, test=adjusted(type=multcomp_adjusted_type))
pq_mppo_miti = summary(mcp_mppo_miti,  
                      test=adjusted(type=multcomp_adjusted_type))$test
mtests_mppo_miti = cbind(pq_mppo_miti$coefficients, pq_mppo_miti$sigma, 
                        pq_mppo_miti$tstat, pq_mppo_miti$pvalues)
error_mppo_miti = attr(pq_mppo_miti$pvalues, "error")
pname_mppo_miti = switch(
  mcp_mppo_miti$alternativ,
  less = paste("Pr(<", ifelse(
    mcp_mppo_miti$dprime_Exp_III_mppo_miti ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(
    mcp_mppo_miti$dprime_Exp_III_mppo_miti == 0, "z", "t"), ")", sep = ""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_mppo_miti$dprime_Exp_III_mppo_miti == 0, "z", "t"), "|)", sep = ""))
colnames(mtests_mppo_miti) = c("Estimate", "Std. Error", ifelse(
  mcp_mppo_miti$df ==0, "z value", "t value"), pname_mppo_miti)

mixmodel.emm = emmeans(mixmodel, ~ MPPO * MITI)
contrast(mixmodel.emm, "eff", by = "MPPO")
contrast(mixmodel.emm, "eff", by = "MITI")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans( mixmodel, c("MPPO","MITI"))

file_name = "interaction_plot_mfpo_miti_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, MPPO ~ MITI)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~MPPO|MITI, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_mfpo_miti_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(mixmodel, ~MPPO|MITI, type="response")

file_name = "pwpp_mfpo_miti_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_mppo_miti)
xtable(emm_df)

################################################################################
# MODELE d' ~ MPPO * TARGET RATE
################################################################################

setwd(exp_path)

dprime_Exp_III_mppo_tra = read.table(
  "psychoac_mppo_targetrate_Exp_III.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_III_mppo_tra = dprime_Exp_III_mppo_tra[which(
  dprime_Exp_III_mppo_tra$Sujet!=8),]

cond_mppo = c("16fpo","32fpo","64fpo")
dprime_Exp_III_mppo_tra$MPPO = factor(
  dprime_Exp_III_mppo_tra$MPPO, levels=cond_mppo)

cond_tra = c("1Hz","2Hz","5Hz")
dprime_Exp_III_mppo_tra$TargetRate = factor(
  dprime_Exp_III_mppo_tra$TargetRate, levels=cond_tra)

qqnorm(dprime_Exp_III_mppo_tra[,"dprime"], main = "QQNORM FOR Dprime")
qqline(dprime_Exp_III_mppo_tra[,"dprime"], main = "QQNORM FOR Dprime")

mixmodel = lme(dprime ~ MPPO * TargetRate, 
               data=dprime_Exp_III_mppo_tra, random=~1|Sujet, method= "ML")

summary(mixmodel)
anova(mixmodel)

dprime_Exp_III_mppo_tra$res = interaction(
  dprime_Exp_III_mppo_tra$MPPO, 
  dprime_Exp_III_mppo_tra$TargetRate)
dprime_Exp_III_mppo_tra$res = factor(dprime_Exp_III_mppo_tra$res)

res_mixmodel = lme(dprime ~ res, data=dprime_Exp_III_mppo_tra, 
                   random=~1|Sujet, method= "ML")

mcp_mppo_tra = glht(res_mixmodel, linfct=mcp(res="Tukey"))
summary(mcp_mppo_tra, test=adjusted(type=multcomp_adjusted_type))
pq_mppo_tra = summary(mcp_mppo_tra,  
                            test=adjusted(type=multcomp_adjusted_type))$test
mtests_mppo_tra = cbind(pq_mppo_tra$coefficients, 
                              pq_mppo_tra$sigma, 
                              pq_mppo_tra$tstat, 
                              pq_mppo_tra$pvalues)
error_mppo_tra = attr(pq_mppo_tra$pvalues, "error")
pname_mppo_tra = switch(
  mcp_mppo_tra$alternativ,
  less = paste("Pr(<", ifelse(
    mcp_mppo_tra$dprime_Exp_III_mppo_tra==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_mppo_tra$dprime_Exp_III_mppo_tra==0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_mppo_tra$dprime_Exp_III_mppo_tra==0, "z", "t"), "|)",sep=""))
colnames(mtests_mppo_tra) = c("Estimate", "Std. Error", ifelse(
  mcp_mppo_tra$df ==0, "z value", "t value"), pname_mppo_tra)

mixmodel.emm = emmeans(mixmodel, ~ MPPO * TargetRate)
contrast(mixmodel.emm, "eff", by = "MPPO")
contrast(mixmodel.emm, "eff", by = "TargetRate")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans( mixmodel, c("MPPO","TargetRate"))

file_name = "interaction_plot_mfpo_targetrate_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, MPPO ~ TargetRate)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~MPPO|TargetRate, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_mfpo_targetrate_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(mixmodel, ~MPPO|TargetRate, type="response")

file_name = "pwpp_mfpo_targetrate_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_mppo_tra)
xtable(emm_df)

################################################################################
# MODELE d' ~ MITI * TARGET RATE
################################################################################

setwd(exp_path)

dprime_Exp_III_miti_tra = read.table(
  "psychoac_miti_targetrate_Exp_III.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_III_miti_tra = dprime_Exp_III_miti_tra[which(
  dprime_Exp_III_miti_tra$Sujet!=8),]

cond_miti = c("2300ms","1100ms","300ms")
dprime_Exp_III_miti_tra$MITI = factor(
  dprime_Exp_III_miti_tra$MITI, levels=cond_miti)

cond_tra = c("1Hz","2Hz","5Hz")
dprime_Exp_III_miti_tra$TargetRate = factor(
  dprime_Exp_III_miti_tra$TargetRate, levels=cond_tra)

qqnorm(dprime_Exp_III_miti_tra[,"dprime"], main = "QQNORM FOR Dprime")
qqline(dprime_Exp_III_miti_tra[,"dprime"], main = "QQNORM FOR Dprime")

mixmodel = lme(dprime ~ MITI * TargetRate, 
               data=dprime_Exp_III_miti_tra, random=~1|Sujet, method= "ML")

summary(mixmodel)
anova(mixmodel)

dprime_Exp_III_miti_tra$res = interaction(
  dprime_Exp_III_miti_tra$MITI, 
  dprime_Exp_III_miti_tra$TargetRate)
dprime_Exp_III_miti_tra$res = factor(dprime_Exp_III_miti_tra$res)

res_mixmodel = lme(dprime ~ res, data=dprime_Exp_III_miti_tra, 
                   random=~1|Sujet, method= "ML")

mcp_miti_tra = glht(res_mixmodel, linfct=mcp(res="Tukey"))
summary(mcp_miti_tra, test=adjusted(type=multcomp_adjusted_type))
pq_miti_tra = summary(mcp_miti_tra,  
                      test=adjusted(type=multcomp_adjusted_type))$test
mtests_miti_tra = cbind(pq_miti_tra$coefficients, 
                        pq_miti_tra$sigma, 
                        pq_miti_tra$tstat, 
                        pq_miti_tra$pvalues)
error_miti_tra = attr(pq_miti_tra$pvalues, "error")
pname_miti_tra = switch(
  mcp_miti_tra$alternativ,
  less = paste("Pr(<", ifelse(
    mcp_miti_tra$dprime_Exp_III_miti_tra==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_miti_tra$dprime_Exp_III_miti_tra==0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_miti_tra$dprime_Exp_III_miti_tra==0, "z", "t"), "|)",sep=""))
colnames(mtests_miti_tra) = c("Estimate", "Std. Error", ifelse(
  mcp_miti_tra$df ==0, "z value", "t value"), pname_miti_tra)

mixmodel.emm = emmeans(mixmodel, ~ MITI * TargetRate)
contrast(mixmodel.emm, "eff", by = "MITI")
contrast(mixmodel.emm, "eff", by = "TargetRate")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans( mixmodel, c("MITI","TargetRate"))

file_name = "interaction_plot_mfpo_targetrate_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, MITI ~ TargetRate)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~MITI|TargetRate, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_mfpo_targetrate_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(mixmodel, ~MITI|TargetRate, type="response")

file_name = "pwpp_mfpo_targetrate_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_miti_tra)
xtable(emm_df)

################################################################################
