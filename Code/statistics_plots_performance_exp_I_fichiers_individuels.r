################################################################################
# SCRIPT FOR PLOTTING AND ANALYZING DPRIME FROM EXPERIMENT I
# THIS PART IS BASED ON MULTIPLE DATAFRAMES WITH ONE FACTOR OR INTERACTION 
# FOR DPRIME FOR EACH SUBJECT
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
article_target = "article_1_PA/PloSOne/maxi_dossier/Exp_I/performance_stats"
slides_path = "/home/link/Documents/thèse_onera/diapos_phd_thesis/"
exp_path = "/home/link/Documents/thèse_onera/experimentation/"
slides_path_acc = "images/EEG/info_content"
slides_path_final = str_c(slides_path, slides_path_acc, sep="")

multcomp_adjusted_type = "fdr"
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")

################################################################################
################################################################################
################################################################################
# THIS PART IS BASED ON SEVERAL DATAFRAMES OF EACH FACTOR AND INTERACTIONS FOR 
# DPRIME FOR EACH SUBJECT
################################################################################
################################################################################
################################################################################

################################################################################
# MODELE d' ~ MPPO
################################################################################

setwd(exp_path)

dprime_Exp_I_mppo = read.table(
  "psychoac_mppo_Exp_I.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I_mppo = dprime_Exp_I_mppo[which(dprime_Exp_I_mppo$Sujet != 2 & 
                                              dprime_Exp_I_mppo$Sujet != 6 & 
                                              dprime_Exp_I_mppo$Sujet != 9),]

cond_mppo = c("4fpo","16fpo","64fpo")
dprime_Exp_I_mppo$MPPO = factor(dprime_Exp_I_mppo$MPPO,levels=cond_mppo)

qqnorm(dprime_Exp_I_mppo[,"dprime"], main = "QQNORM FOR Dprime MPPO")
qqline(dprime_Exp_I_mppo[,"dprime"], main = "QQNORM FOR Dprime MPPO")

mixmodel = lme(dprime~MPPO, data=dprime_Exp_I_mppo, 
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

dprime_Exp_I_mtd = read.table(
  "psychoac_mtd_Exp_I.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I_mtd = dprime_Exp_I_mtd[which(dprime_Exp_I_mtd$Sujet != 2 & 
                                            dprime_Exp_I_mtd$Sujet != 6 & 
                                            dprime_Exp_I_mtd$Sujet != 9),]

cond_mtd = c("20ms","60ms","100ms")
dprime_Exp_I_mtd$MTD = factor(dprime_Exp_I_mtd$MTD, levels=cond_mtd)

qqnorm(dprime_Exp_I_mtd[,"dprime"], main = "QQNORM FOR Dprime MTD")
qqline(dprime_Exp_I_mtd[,"dprime"], main = "QQNORM FOR Dprime MTD")

mixmodel = lme(dprime~MTD, data=dprime_Exp_I_mtd, random=~1|Sujet, method="ML")
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
# MODELE d' ~ TTD
################################################################################

setwd(exp_path)

dprime_Exp_I_ttd = read.table(
  "psychoac_ttd_Exp_I.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I_ttd = dprime_Exp_I_ttd[which(dprime_Exp_I_ttd$Sujet != 2 & 
                                            dprime_Exp_I_ttd$Sujet != 6 & 
                                            dprime_Exp_I_ttd$Sujet != 9),]

cond_ttd = c("20ms","60ms","100ms")
dprime_Exp_I_ttd$TTD = factor(dprime_Exp_I_ttd$TTD, levels=cond_ttd)

qqnorm(dprime_Exp_I_ttd[,"dprime"], main = "QQNORM FOR Dprime TTD")
qqline(dprime_Exp_I_ttd[,"dprime"], main = "QQNORM FOR Dprime TTD")

mixmodel = lme(dprime~TTD, data=dprime_Exp_I_ttd, random=~1|Sujet, method="ML")
summary(mixmodel)
anova(mixmodel)

mcp_ttd = glht(mixmodel, linfct=mcp(TTD="Tukey"))
summary(mcp_ttd,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_ttd)$test
mtests_ttd = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_ttd$alternativ,
  less = paste("Pr(<", ifelse(mcp_ttd$df ==0, "z", "t"), ")", sep = ""), 
  greater = paste("Pr(>", ifelse(mcp_ttd$df == 0, "z", "t"), ")", sep = ""), 
  two.sided = paste("Pr(>|", ifelse(mcp_ttd$df == 0, "z", "t"), "|)", sep = ""))
colnames(mtests_ttd) = c("Estimate", "Std. Error", 
                         ifelse(mcp_mtd$df ==0, "z value", "t value"), pname)
report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_ttd)

################################################################################
# MODELE d' ~ UNCERTAINTY
################################################################################

setwd(exp_path)

dprime_Exp_I_uncertainty = read.table(
  "psychoac_uncertainty_Exp_I.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I_uncertainty = 
  dprime_Exp_I_uncertainty[which(dprime_Exp_I_uncertainty$Sujet != 2 & 
                                   dprime_Exp_I_uncertainty$Sujet != 6 & 
                                   dprime_Exp_I_uncertainty$Sujet != 9),]

dprime_Exp_I_uncertainty$Uncertainty = as.factor(
  dprime_Exp_I_uncertainty$Uncertainty)

qqnorm(dprime_Exp_I_uncertainty[,"dprime"], 
       main = "QQNORM FOR Dprime Uncertainty")
qqline(dprime_Exp_I_uncertainty[,"dprime"], 
       main = "QQNORM FOR Dprime Uncertainty")

mixmodel = lme(dprime~Uncertainty, data=dprime_Exp_I_uncertainty, 
               random=~1|Sujet, method="ML")
summary(mixmodel)
anova(mixmodel)

mcp_uncertainty = glht(mixmodel, linfct=mcp(Uncertainty="Tukey"))
summary(mcp_uncertainty,test=adjusted(type=multcomp_adjusted_type))
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

dprime_Exp_I_similarity = read.table(
  "psychoac_similarity_Exp_I.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I_similarity = 
  dprime_Exp_I_similarity[which(dprime_Exp_I_similarity$Sujet != 2 & 
                                  dprime_Exp_I_similarity$Sujet != 6 & 
                                  dprime_Exp_I_similarity$Sujet != 9),]

cond_similarity = c("-80","-40","0","40","80")
dprime_Exp_I_similarity$Similarity = factor(
  dprime_Exp_I_similarity$Similarity, levels=cond_similarity)

qqnorm(dprime_Exp_I_similarity[,"dprime"], 
       main = "QQNORM FOR Dprime Similarity")
qqline(dprime_Exp_I_similarity[,"dprime"], 
       main = "QQNORM FOR Dprime Similarity")

mixmodel = lme(dprime~Similarity, data=dprime_Exp_I_similarity, 
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
# MODELE d' ~ MPPO * MTD
################################################################################

setwd(exp_path)

dprime_Exp_I_mppo_mtd = read.table(
  "psychoac_mppo_mtd_Exp_I.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I_mppo_mtd = 
  dprime_Exp_I_mppo_mtd[which(dprime_Exp_I_mppo_mtd$Sujet != 2 & 
                                dprime_Exp_I_mppo_mtd$Sujet != 6 & 
                                dprime_Exp_I_mppo_mtd$Sujet != 9),]

cond_mppo = c("4fpo","16fpo","64fpo")
dprime_Exp_I_mppo_mtd$MPPO = factor(dprime_Exp_I_mppo_mtd$MPPO,levels=cond_mppo)

cond_mtd = c("20ms","60ms","100ms")
dprime_Exp_I_mppo_mtd$MTD = factor(dprime_Exp_I_mppo_mtd$MTD, levels=cond_mtd)

qqnorm(dprime_Exp_I_mppo_mtd[,"dprime"], main = "QQNORM FOR Dprime")
qqline(dprime_Exp_I_mppo_mtd[,"dprime"], main = "QQNORM FOR Dprime")

mixmodel = lme(dprime ~ MPPO * MTD, data=dprime_Exp_I_mppo_mtd,
               random=~1|Sujet, method= "ML")

summary(mixmodel)
anova(mixmodel)

dprime_Exp_I_mppo_mtd$res = interaction(dprime_Exp_I_mppo_mtd$MPPO, 
                                        dprime_Exp_I_mppo_mtd$MTD)
dprime_Exp_I_mppo_mtd$res = factor(dprime_Exp_I_mppo_mtd$res)

res_mixmodel = lme(dprime ~ res, data=dprime_Exp_I_mppo_mtd, random=~1|Sujet, 
                   method= "ML")

mcp_mppo_mtd = glht(res_mixmodel, linfct=mcp(res="Tukey"))
summary(mcp_mppo_mtd, test=adjusted(type=multcomp_adjusted_type))
pq_mppo_mtd = summary(mcp_mppo_mtd,  
                      test=adjusted(type=multcomp_adjusted_type))$test
mtests_mppo_mtd = cbind(pq_mppo_mtd$coefficients, pq_mppo_mtd$sigma, 
                        pq_mppo_mtd$tstat, pq_mppo_mtd$pvalues)
error_mppo_mtd = attr(pq_mppo_mtd$pvalues, "error")
pname_mppo_mtd = switch(
  mcp_mppo_mtd$alternativ,
  less = paste("Pr(<", ifelse(
    mcp_mppo_mtd$dprime_Exp_I_mppo_mtd ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(
    mcp_mppo_mtd$dprime_Exp_I_mppo_mtd == 0, "z", "t"), ")", sep = ""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_mppo_mtd$dprime_Exp_I_mppo_mtd == 0, "z", "t"), "|)", sep = ""))
colnames(mtests_mppo_mtd) = c("Estimate", "Std. Error", ifelse(
  mcp_mppo_mtd$df ==0, "z value", "t value"), pname_mppo_mtd)

mixmodel.emm = emmeans(mixmodel, ~ MPPO * MTD)
contrast(mixmodel.emm, "eff", by = "MPPO")
contrast(mixmodel.emm, "eff", by = "MTD")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans( mixmodel, c("MPPO","MTD"))

file_name = "interaction_plot_mfpo_mtd_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, MPPO ~ MTD)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~MPPO|MTD, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_mfpo_mtd_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(mixmodel, ~MPPO|MTD, type="response")

file_name = "pwpp_mfpo_mtd_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

report(mixmodel)
xtable(anova(mixmodel))
xtable(mtests_mppo_mtd)
xtable(emm_df)

################################################################################
# MODELE d' ~ UNCERTAINTY * SIMILARITY
################################################################################

setwd(exp_path)

dprime_Exp_I_uncertainty_similarity = read.table(
  "psychoac_uncertainty_similarity_Exp_I.csv", 
  header=TRUE, sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I_uncertainty_similarity = 
  dprime_Exp_I_uncertainty_similarity[which(
    dprime_Exp_I_uncertainty_similarity$Sujet != 2 & 
      dprime_Exp_I_uncertainty_similarity$Sujet != 6 & 
      dprime_Exp_I_uncertainty_similarity$Sujet != 9),]

cond_uncertainty = c("29","115","463")
dprime_Exp_I_uncertainty_similarity$Uncertainty = factor(
  dprime_Exp_I_uncertainty_similarity$Uncertainty,levels=cond_uncertainty)

cond_similarity = c("-80","-40","0","40","80")
dprime_Exp_I_uncertainty_similarity$Similarity = factor(
  dprime_Exp_I_uncertainty_similarity$Similarity, levels=cond_similarity)

qqnorm(dprime_Exp_I_uncertainty_similarity[,"dprime"], main="QQNORM FOR Dprime")
qqline(dprime_Exp_I_uncertainty_similarity[,"dprime"], main="QQNORM FOR Dprime")

mixmodel = lme(dprime ~ Uncertainty * Similarity, 
               data=dprime_Exp_I_uncertainty_similarity, 
               random=~1|Sujet, method= "ML")

summary(mixmodel)
anova(mixmodel)

dprime_Exp_I_uncertainty_similarity$res = interaction(
  dprime_Exp_I_uncertainty_similarity$Uncertainty, 
  dprime_Exp_I_uncertainty_similarity$Similarity)
dprime_Exp_I_uncertainty_similarity$res = factor(
  dprime_Exp_I_uncertainty_similarity$res)

res_mixmodel = lme(dprime ~ res, data=dprime_Exp_I_uncertainty_similarity, 
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
    mcp_uncertainty_similarity$dprime_Exp_I_uncertainty_similarity == 0, 
    "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(
    mcp_uncertainty_similarity$dprime_Exp_I_uncertainty_similarity == 0, 
    "z", "t"), ")", sep = ""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_uncertainty_similarity$dprime_Exp_I_uncertainty_similarity == 0, 
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

file_name = "interaction_plot_uncertainty_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, Uncertainty ~ Similarity)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~Uncertainty|Similarity, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_uncertainty_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(mixmodel, ~Uncertainty|Similarity, type="response")

file_name = "pwpp_uncertainty_similarity_exp_I.pdf"
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
