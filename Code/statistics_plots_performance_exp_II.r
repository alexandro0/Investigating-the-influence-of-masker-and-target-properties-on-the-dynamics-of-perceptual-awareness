################################################################################
# SCRIPT FOR PLOTTING AND ANALYZING DPRIME FROM EXPERIMENT II
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
article_target = "article_1_PA/PloSOne/maxi_dossier/Exp_II/performance_stats"
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

dprime_Exp_II_mppo = read.table(
  "psychoac_mppo_Exp_II.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_II_mppo = dprime_Exp_II_mppo[which(
  dprime_Exp_II_mppo$Sujet != 1 & dprime_Exp_II_mppo$Sujet != 2 & 
    dprime_Exp_II_mppo$Sujet != 6 & dprime_Exp_II_mppo$Sujet != 8),]

cond_mppo = c("4fpo","16fpo","64fpo")
dprime_Exp_II_mppo$MPPO = factor(dprime_Exp_II_mppo$MPPO,levels=cond_mppo)

qqnorm(dprime_Exp_II_mppo[,"dprime"], main = "QQNORM FOR Dprime MPPO")
qqline(dprime_Exp_II_mppo[,"dprime"], main = "QQNORM FOR Dprime MPPO")

mixmodel = lme(dprime~MPPO, data=dprime_Exp_II_mppo, 
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
# MODELE d' ~ MTD
################################################################################

setwd(exp_path)

dprime_Exp_II_mtd = read.table(
  "psychoac_mtd_Exp_II.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_II_mtd = dprime_Exp_II_mtd[which(
  dprime_Exp_II_mtd$Sujet != 1 & dprime_Exp_II_mtd$Sujet != 2 & 
    dprime_Exp_II_mtd$Sujet != 6 & dprime_Exp_II_mtd$Sujet != 8),]

cond_mtd = c("20ms","60ms","100ms")
dprime_Exp_II_mtd$MTD = factor(dprime_Exp_II_mtd$MTD, levels=cond_mtd)

qqnorm(dprime_Exp_II_mtd[,"dprime"], main = "QQNORM FOR Dprime MTD")
qqline(dprime_Exp_II_mtd[,"dprime"], main = "QQNORM FOR Dprime MTD")

mixmodel = lme(dprime~MTD, data=dprime_Exp_II_mtd, random=~1|Sujet, method="ML")
summary(mixmodel)
anova(mixmodel)

mcp_mtd = glht(mixmodel, linfct=mcp(MTD="Tukey"))
summary(mcp_mtd,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_mtd)$test
mtests_mtd = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_mtd$alternativ,
  less = paste("Pr(<", ifelse(mcp_mtd$df ==0, "z", "t"), ")", sep = ""), 
  greater = paste("Pr(>", ifelse(mcp_mtd$df == 0, "z", "t"), ")", sep = ""), 
  two.sided = paste("Pr(>|", ifelse(mcp_mtd$df == 0, "z", "t"), "|)", sep = ""))                                                                   
colnames(mtests_mtd) = c("Estimate", "Std. Error", 
                         ifelse(mcp_mtd$df ==0, "z value", "t value"), pname)
report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_mtd)

################################################################################
# MODELE d' ~ TARGET RATE
################################################################################

setwd(exp_path)

dprime_Exp_II_tra = read.table(
  "psychoac_targetrate_Exp_II.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_II_tra = dprime_Exp_II_tra[which(
  dprime_Exp_II_tra$Sujet != 1 & dprime_Exp_II_tra$Sujet != 2 & 
    dprime_Exp_II_tra$Sujet != 6 & dprime_Exp_II_tra$Sujet != 8),]

cond_tra = c("5Hz","10Hz","20Hz")
dprime_Exp_II_tra$TargetRate = factor(dprime_Exp_II_tra$TargetRate, 
                                      levels=cond_tra)

qqnorm(dprime_Exp_II_tra[,"dprime"], main = "QQNORM FOR Dprime TargetRate")
qqline(dprime_Exp_II_tra[,"dprime"], main = "QQNORM FOR Dprime TargetRate")

mixmodel = lme(dprime~TargetRate, data=dprime_Exp_II_tra, 
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

dprime_Exp_II_uncertainty = read.table(
  "psychoac_uncertainty_Exp_II.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_II_uncertainty = dprime_Exp_II_uncertainty[which(
  dprime_Exp_II_uncertainty$Sujet != 1 & 
    dprime_Exp_II_uncertainty$Sujet != 2 & 
    dprime_Exp_II_uncertainty$Sujet != 6 & 
    dprime_Exp_II_uncertainty$Sujet != 8),]

cond_uncertainty = c("29","115","463")
dprime_Exp_II_uncertainty$Uncertainty = factor(
  dprime_Exp_II_uncertainty$Uncertainty,levels=cond_uncertainty)

qqnorm(dprime_Exp_II_uncertainty[,"dprime"], 
       main = "QQNORM FOR Dprime Uncertainty")
qqline(dprime_Exp_II_uncertainty[,"dprime"], 
       main = "QQNORM FOR Dprime Uncertainty")

mixmodel = lme(dprime~Uncertainty, data=dprime_Exp_II_uncertainty, 
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
# MODELE d' ~ SIMILARITY
################################################################################

setwd(exp_path)

dprime_Exp_II_similarity = read.table(
  "psychoac_similarity_Exp_II.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_II_similarity = dprime_Exp_II_similarity[which(
  dprime_Exp_II_similarity$Sujet != 1 & 
    dprime_Exp_II_similarity$Sujet != 2 & 
    dprime_Exp_II_similarity$Sujet != 6 & 
    dprime_Exp_II_similarity$Sujet != 8),]

cond_similarity = c("0","40","80")
dprime_Exp_II_similarity$Similarity = factor(
  dprime_Exp_II_similarity$Similarity, levels=cond_similarity)

qqnorm(dprime_Exp_II_similarity[,"dprime"], 
       main = "QQNORM FOR Dprime Similarity")
qqline(dprime_Exp_II_similarity[,"dprime"], 
       main = "QQNORM FOR Dprime Similarity")

mixmodel = lme(dprime~Similarity, data=dprime_Exp_II_similarity, 
               random=~1|Sujet, method="ML")
summary(mixmodel)
anova(mixmodel)

mcp_similarity = glht(mixmodel, linfct=mcp(Similarity="Tukey"))
summary(mcp_similarity,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_similarity)$test
mtests_similarity = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_similarity$alternativ, 
  less = paste("Pr(<", ifelse(mcp_similarity$df ==0,"z","t"),")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_similarity$df ==0,"z","t"),")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_similarity$df ==0,"z","t"),"|)",sep=""))
colnames(mtests_similarity) = c("Estimate", "Std. Error", 
                                ifelse(mcp_similarity$df ==0, 
                                       "z value", "t value"), pname)
report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_similarity)

################################################################################
# MODELE d' ~ SIMILARITY * TARGET RATE
################################################################################

setwd(exp_path)

dprime_Exp_II_similarity_tra = read.table(
  "psychoac_similarity_targetrate_Exp_II.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_II_similarity_tra = dprime_Exp_II_similarity_tra[which(
  dprime_Exp_II_similarity_tra$Sujet != 1 & 
    dprime_Exp_II_similarity_tra$Sujet != 2 & 
    dprime_Exp_II_similarity_tra$Sujet != 6 & 
    dprime_Exp_II_similarity_tra$Sujet != 8),]

cond_similarity = c("0","40","80")
dprime_Exp_II_similarity_tra$Similarity = factor(
  dprime_Exp_II_similarity_tra$Similarity,levels=cond_similarity)

cond_tra = c("5Hz","10Hz","20Hz")
dprime_Exp_II_similarity_tra$TargetRate = factor(
  dprime_Exp_II_similarity_tra$TargetRate, levels=cond_tra)

qqnorm(dprime_Exp_II_similarity_tra[,"dprime"], main = "QQNORM FOR Dprime")
qqline(dprime_Exp_II_similarity_tra[,"dprime"], main = "QQNORM FOR Dprime")

mixmodel = lme(dprime ~ Similarity * TargetRate, 
               data=dprime_Exp_II_similarity_tra, random=~1|Sujet, method= "ML")

summary(mixmodel)
anova(mixmodel)

dprime_Exp_II_similarity_tra$res = interaction(
  dprime_Exp_II_similarity_tra$Similarity, 
  dprime_Exp_II_similarity_tra$TargetRate)
dprime_Exp_II_similarity_tra$res = factor(dprime_Exp_II_similarity_tra$res)

res_mixmodel = lme(dprime ~ res, data=dprime_Exp_II_similarity_tra, 
                   random=~1|Sujet, method= "ML")

mcp_similarity_tra = glht(res_mixmodel, linfct=mcp(res="Tukey"))
summary(mcp_similarity_tra, test=adjusted(type=multcomp_adjusted_type))
pq_similarity_tra = summary(mcp_similarity_tra,  
                      test=adjusted(type=multcomp_adjusted_type))$test
mtests_similarity_tra = cbind(pq_similarity_tra$coefficients, 
                              pq_similarity_tra$sigma, 
                              pq_similarity_tra$tstat, 
                              pq_similarity_tra$pvalues)
error_similarity_tra = attr(pq_similarity_tra$pvalues, "error")
pname_similarity_tra = switch(
  mcp_similarity_tra$alternativ,
  less = paste("Pr(<", ifelse(
    mcp_similarity_tra$dprime_Exp_II_similarity_tra==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_similarity_tra$dprime_Exp_II_similarity_tra==0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_similarity_tra$dprime_Exp_II_similarity_tra==0, "z", "t"), "|)",sep=""))
colnames(mtests_similarity_tra) = c("Estimate", "Std. Error", ifelse(
  mcp_similarity_tra$df ==0, "z value", "t value"), pname_similarity_tra)

mixmodel.emm = emmeans(mixmodel, ~ Similarity * TargetRate)
contrast(mixmodel.emm, "eff", by = "Similarity")
contrast(mixmodel.emm, "eff", by = "TargetRate")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans( mixmodel, c("Similarity","TargetRate"))

file_name = "interaction_plot_similarity_targetrate_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, Similarity ~ TargetRate)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~Similarity|TargetRate, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_similarity_targetrate_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(mixmodel, ~Similarity|TargetRate, type="response")

file_name = "pwpp_similarity_targetrate_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_similarity_tra)
xtable(emm_df)

################################################################################
# MODELE d' ~ UNCERTAINTY * SIMILARITY
################################################################################

setwd(exp_path)

dprime_Exp_II_uncertainty_similarity = read.table(
  "psychoac_uncertainty_similarity_Exp_II.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_II_uncertainty_similarity = 
  dprime_Exp_II_uncertainty_similarity[which(
    dprime_Exp_II_uncertainty_similarity$Sujet != 1 & 
      dprime_Exp_II_uncertainty_similarity$Sujet != 2 & 
      dprime_Exp_II_uncertainty_similarity$Sujet != 6 & 
      dprime_Exp_II_uncertainty_similarity$Sujet != 8),]

cond_uncertainty = c("29","115","463")
dprime_Exp_II_uncertainty_similarity$Uncertainty = factor(
  dprime_Exp_II_uncertainty_similarity$Uncertainty,levels=cond_uncertainty)

cond_similarity = c("0","40","80")
dprime_Exp_II_uncertainty_similarity$Similarity = factor(
  dprime_Exp_II_uncertainty_similarity$Similarity, levels=cond_similarity)

qqnorm(dprime_Exp_II_uncertainty_similarity[,"dprime"],main="QQNORM FOR Dprime")
qqline(dprime_Exp_II_uncertainty_similarity[,"dprime"],main="QQNORM FOR Dprime")

mixmodel = lme(dprime ~ Uncertainty * Similarity, 
               data=dprime_Exp_II_uncertainty_similarity, 
               random=~1|Sujet, method= "ML")

summary(mixmodel)
anova(mixmodel)

dprime_Exp_II_uncertainty_similarity$res = interaction(
  dprime_Exp_II_uncertainty_similarity$Uncertainty, 
  dprime_Exp_II_uncertainty_similarity$Similarity)
dprime_Exp_II_uncertainty_similarity$res = factor(
  dprime_Exp_II_uncertainty_similarity$res)

res_mixmodel = lme(dprime ~ res, data=dprime_Exp_II_uncertainty_similarity, 
                   random=~1|Sujet, method= "ML")

mcp_uncertainty_similarity = glht(res_mixmodel, linfct=mcp(res="Tukey"))
summary(mcp_uncertainty_similarity, test=adjusted(type=multcomp_adjusted_type))
pq_uncertainty_similarity = summary(mcp_uncertainty_similarity,  
                                    test=adjusted(
                                      type=multcomp_adjusted_type))$test
mtests_uncertainty_similarity = cbind(pq_uncertainty_similarity$coefficients, 
                                      pq_uncertainty_similarity$sigma, 
                                      pq_uncertainty_similarity$tstat, 
                                      pq_uncertainty_similarity$pvalues)
error_uncertainty_similarity = attr(pq_uncertainty_similarity$pvalues, "error")
pname_uncertainty_similarity = switch(
  mcp_uncertainty_similarity$alternativ,
  less = paste("Pr(<", ifelse(
    mcp_uncertainty_similarity$dprime_Exp_II_uncertainty_similarity == 0, 
    "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(
    mcp_uncertainty_similarity$dprime_Exp_II_uncertainty_similarity == 0, 
    "z", "t"), ")", sep = ""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_uncertainty_similarity$dprime_Exp_II_uncertainty_similarity == 0, 
    "z", "t"), "|)", sep = ""))
colnames(mtests_uncertainty_similarity) = c("Estimate", "Std. Error", ifelse(
  mcp_uncertainty_similarity$df ==0, "z value", "t value"), 
  pname_uncertainty_similarity)

mixmodel.emm = emmeans(mixmodel, ~ Uncertainty * Similarity)
contrast(mixmodel.emm, "eff", by = "Uncertainty")
contrast(mixmodel.emm, "eff", by = "Similarity")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans( mixmodel, c("Uncertainty","Similarity"))

file_name = "interaction_plot_uncertainty_similarity_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, Uncertainty ~ Similarity)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~Uncertainty|Similarity, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_uncertainty_similarity_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(mixmodel, ~Uncertainty|Similarity, type="response")

file_name = "pwpp_uncertainty_similarity_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_uncertainty_similarity)
xtable(emm_df)

################################################################################
# MODELE d' ~ UNCERTAINTY * TARGET RATE
################################################################################

setwd(exp_path)

dprime_Exp_II_uncertainty_tra = read.table(
  "psychoac_uncertainty_targetrate_Exp_II.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_II_uncertainty_tra = 
  dprime_Exp_II_uncertainty_tra[which(
    dprime_Exp_II_uncertainty_tra$Sujet != 1 & 
      dprime_Exp_II_uncertainty_tra$Sujet != 2 & 
      dprime_Exp_II_uncertainty_tra$Sujet != 6 & 
      dprime_Exp_II_uncertainty_tra$Sujet != 8),]

cond_uncertainty = c("29","115","463")
dprime_Exp_II_uncertainty_tra$Uncertainty = factor(
  dprime_Exp_II_uncertainty_tra$Uncertainty,levels=cond_uncertainty)

cond_tra = c("5Hz","10Hz","20Hz")
dprime_Exp_II_uncertainty_tra$TargetRate = factor(
  dprime_Exp_II_uncertainty_tra$TargetRate, levels=cond_tra)

qqnorm(dprime_Exp_II_uncertainty_tra[,"dprime"],main="QQNORM FOR Dprime")
qqline(dprime_Exp_II_uncertainty_tra[,"dprime"],main="QQNORM FOR Dprime")

mixmodel = lme(dprime ~ Uncertainty * TargetRate, 
               data=dprime_Exp_II_uncertainty_tra, 
               random=~1|Sujet, method= "ML")

summary(mixmodel)
anova(mixmodel)

dprime_Exp_II_uncertainty_tra$res = interaction(
  dprime_Exp_II_uncertainty_tra$Uncertainty, 
  dprime_Exp_II_uncertainty_tra$TargetRate)
dprime_Exp_II_uncertainty_tra$res = factor(
  dprime_Exp_II_uncertainty_tra$res)

res_mixmodel = lme(dprime ~ res, data=dprime_Exp_II_uncertainty_tra, 
                   random=~1|Sujet, method= "ML")

mcp_uncertainty_tra = glht(res_mixmodel, linfct=mcp(res="Tukey"))
summary(mcp_uncertainty_tra, test=adjusted(type=multcomp_adjusted_type))
pq_uncertainty_tra = summary(mcp_uncertainty_tra,  
                                    test=adjusted(
                                      type=multcomp_adjusted_type))$test
mtests_uncertainty_tra = cbind(pq_uncertainty_tra$coefficients, 
                                      pq_uncertainty_tra$sigma, 
                                      pq_uncertainty_tra$tstat, 
                                      pq_uncertainty_tra$pvalues)
error_uncertainty_tra = attr(pq_uncertainty_tra$pvalues, "error")
pname_uncertainty_tra = switch(
  mcp_uncertainty_tra$alternativ,
  less = paste("Pr(<", ifelse(
    mcp_uncertainty_tra$dprime_Exp_II_uncertainty_tra == 0, 
    "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(
    mcp_uncertainty_tra$dprime_Exp_II_uncertainty_tra == 0, 
    "z", "t"), ")", sep = ""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_uncertainty_tra$dprime_Exp_II_uncertainty_tra == 0, 
    "z", "t"), "|)", sep = ""))
colnames(mtests_uncertainty_tra) = c("Estimate", "Std. Error", ifelse(
  mcp_uncertainty_tra$df ==0, "z value", "t value"), 
  pname_uncertainty_tra)

mixmodel.emm = emmeans(mixmodel, ~ Uncertainty * TargetRate)
contrast(mixmodel.emm, "eff", by = "Uncertainty")
contrast(mixmodel.emm, "eff", by = "TargetRate")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans( mixmodel, c("Uncertainty","TargetRate"))

file_name = "interaction_plot_uncertainty_targetrate_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, Uncertainty ~ TargetRate)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~Uncertainty|TargetRate, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_uncertainty_targetrate_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(mixmodel, ~Uncertainty|TargetRate, type="response")

file_name = "pwpp_uncertainty_targetrate_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_uncertainty_tra)
xtable(emm_df)

################################################################################
