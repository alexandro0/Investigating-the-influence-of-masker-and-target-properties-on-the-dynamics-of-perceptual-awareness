################################################################################
# SCRIPT FOR PLOTTING AND ANALYZING DPRIME FROM EXPERIMENT I
# THIS PART IS BASED ON DATAFRAME WITH ALL FACTORS AND INTERACTIONS FOR DPRIME 
# FOR EACH SUBJECT
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
# MODELE DPRIME ~ MFPO * MTD * TTD
################################################################################

setwd(exp_path)

dprime_Exp_I = read.table("psychoac_all_factors_Exp_I.csv", header=TRUE, 
                               sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I = dprime_Exp_I[which(dprime_Exp_I$Sujet != 2 & 
                                    dprime_Exp_I$Sujet != 6 & 
                                    dprime_Exp_I$Sujet != 9),]

cond_mppo = c("4fpo","16fpo","64fpo")
dprime_Exp_I$MFPO = factor(dprime_Exp_I$MFPO,levels=cond_mppo)

cond_mtd = c("20ms","60ms","100ms")
dprime_Exp_I$MTD = factor(dprime_Exp_I$MTD, levels=cond_mtd)

cond_ttd = c("20ms","60ms","100ms")
dprime_Exp_I$TTD = factor(dprime_Exp_I$TTD, levels=cond_ttd)

mixmodel = lme(dprime ~ MFPO * MTD * TTD, data=dprime_Exp_I, random=~1|Sujet, 
               method="ML")

summary(mixmodel)
anova(mixmodel)

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# MFPO
mcp_mppo = glht(mixmodel, linfct=mcp(MFPO="Tukey", interaction_average = TRUE,
                                  covariate_average = TRUE))
summary(mcp_mppo, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_mppo)$test
mtests_mppo_exp_I=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_mppo$alternativ,
  less=paste("Pr(<", ifelse(mcp_mppo$df ==0, "z", "t"), ")", sep = ""),
  greater=paste("Pr(>", ifelse(mcp_mppo$df==0, "z", "t"), ")", sep = ""),
  two.sided=paste("Pr(>|", ifelse(mcp_mppo$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_mppo_exp_I) = c("Estimate", "Std. Error",
                                ifelse(mcp_mppo$df ==0, "z value", "t value"),
                                pname)

# MTD
mcp_mtd = glht(mixmodel, linfct=mcp(MTD="Tukey", interaction_average = TRUE,
                                 covariate_average = TRUE))
summary(mcp_mtd, test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_mtd)$test
mtests_mtd_exp_I = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_mtd$alternativ,
  less = paste("Pr(<", ifelse(mcp_mtd$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_mtd$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_mtd$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_mtd_exp_I) = c("Estimate", "Std. Error",
                               ifelse(mcp_mtd$df==0, "z value", "t value"),
                               pname)

# TTD
mcp_ttd = glht(mixmodel, linfct=mcp(TTD="Tukey", interaction_average = TRUE,
                                 covariate_average = TRUE))
summary(mcp_ttd,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_ttd)$test
mtests_ttd_exp_I = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_ttd$alternativ,
  less = paste("Pr(<", ifelse(mcp_ttd$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_ttd$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_ttd$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_ttd_exp_I) = c("Estimate", "Std. Error",
                               ifelse(mcp_ttd$df==0, "z value", "t value"),
                               pname)

# MFPO * MTD
dprime_Exp_I$res = interaction(dprime_Exp_I$MFPO, dprime_Exp_I$MTD)
dprime_Exp_I$res = factor(dprime_Exp_I$res)
model_mppo_mtd = lme(dprime ~ res, data=dprime_Exp_I, random=~1|Sujet, 
                     method="ML")
mcp_mppo_mtd = glht(model_mppo_mtd, linfct=mcp(res="Tukey"))
summary(mcp_mppo_mtd, test=adjusted(type="bonferroni"))
pq_mppo_mtd = summary(mcp_mppo_mtd, test=adjusted(type="bonferroni"))$test
mtests_mppo_mtd_exp_I = cbind(pq_mppo_mtd$coefficients, pq_mppo_mtd$sigma, 
                              pq_mppo_mtd$tstat, pq_mppo_mtd$pvalues)
error_mppo_mtd = attr(pq_mppo_mtd$pvalues, "error")
pname_mppo_mtd = switch(
  mcp_mppo_mtd$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_mtd$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_mtd$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_mtd$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_mtd_exp_I) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_mtd$df==0,"z value","t value"), pname_mppo_mtd)
# mtests_mppo_mtd_exp_I = as.data.frame(mtests_mppo_mtd_exp_I)
# mtests_mppo_mtd_exp_I = mtests_mppo_mtd_exp_I[
#   which(mtests_mppo_mtd_exp_I$`Pr(>|z|)`<0.05),]

# MFPO * TTD
dprime_Exp_I$res = interaction(dprime_Exp_I$MFPO, dprime_Exp_I$TTD)
dprime_Exp_I$res = factor(dprime_Exp_I$res)
model_mppo_ttd = lme(dprime ~ res, data=dprime_Exp_I, random=~1|Sujet, 
                     method="ML")
mcp_mppo_ttd = glht(model_mppo_ttd, linfct=mcp(res="Tukey"))
summary(mcp_mppo_ttd, test=adjusted(type="bonferroni"))
pq_mppo_ttd = summary(mcp_mppo_ttd, test=adjusted(type="bonferroni"))$test
mtests_mppo_ttd_exp_I = cbind(pq_mppo_ttd$coefficients, pq_mppo_ttd$sigma, 
                              pq_mppo_ttd$tstat, pq_mppo_ttd$pvalues)
error_mppo_ttd = attr(pq_mppo_ttd$pvalues, "error")
pname_mppo_ttd = switch(
  mcp_mppo_ttd$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_ttd$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_ttd$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_ttd$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_ttd_exp_I) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_ttd$df==0,"z value","t value"), pname_mppo_ttd)
# mtests_mppo_ttd_exp_I = as.data.frame(mtests_mppo_ttd_exp_I)
# mtests_mppo_ttd_exp_I = mtests_mppo_ttd_exp_I[
#   which(mtests_mppo_ttd_exp_I$`Pr(>|z|)`<0.05),]

# MTD * TTD
dprime_Exp_I$res = interaction(dprime_Exp_I$MTD, dprime_Exp_I$TTD)
dprime_Exp_I$res = factor(dprime_Exp_I$res)
model_mtd_ttd = lme(dprime ~ res, data=dprime_Exp_I, random=~1|Sujet, 
                     method="ML")
mcp_mtd_ttd = glht(model_mtd_ttd, linfct=mcp(res="Tukey"))
summary(mcp_mtd_ttd, test=adjusted(type="bonferroni"))
pq_mtd_ttd = summary(mcp_mtd_ttd, test=adjusted(type="bonferroni"))$test
mtests_mtd_ttd_exp_I = cbind(pq_mtd_ttd$coefficients, pq_mtd_ttd$sigma, 
                             pq_mtd_ttd$tstat, pq_mtd_ttd$pvalues)
error_mtd_ttd = attr(pq_mtd_ttd$pvalues, "error")
pname_mtd_ttd = switch(
  mcp_mtd_ttd$alternativ,
  less = paste("Pr(<", ifelse(mcp_mtd_ttd$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mtd_ttd$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mtd_ttd$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mtd_ttd_exp_I) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mtd_ttd$df==0,"z value","t value"), pname_mtd_ttd)
# mtests_mtd_ttd_exp_I = as.data.frame(mtests_mtd_ttd_exp_I)
# mtests_mtd_ttd_exp_I = mtests_mtd_ttd_exp_I[
#   which(mtests_mtd_ttd_exp_I$`Pr(>|z|)`<0.05),]

# MFPO * MTD * TTD
dprime_Exp_I$res = interaction(
  dprime_Exp_I$MFPO, dprime_Exp_I$MTD, dprime_Exp_I$TTD)
dprime_Exp_I$res = factor(dprime_Exp_I$res)
model_mppo_mtd_ttd = lme(dprime ~ res, data=dprime_Exp_I, random=~1|Sujet, 
                         method="ML")
mcp_mppo_mtd_ttd = glht(model_mppo_mtd_ttd, linfct=mcp(res="Tukey"))
summary(mcp_mppo_mtd_ttd, test=adjusted(type="bonferroni"))
pq_mppo_mtd_ttd = summary(mcp_mppo_mtd_ttd, 
                          test=adjusted(type="bonferroni"))$test
mtests_mppo_mtd_ttd_exp_I = cbind(
  pq_mppo_mtd_ttd$coefficients, pq_mppo_mtd_ttd$sigma, 
  pq_mppo_mtd_ttd$tstat, pq_mppo_mtd_ttd$pvalues)
error_mppo_mtd_ttd = attr(pq_mppo_mtd_ttd$pvalues, "error")
pname_mppo_mtd_ttd = switch(
  mcp_mppo_mtd_ttd$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_mtd_ttd$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_mppo_mtd_ttd$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_mppo_mtd_ttd$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_mtd_ttd_exp_I) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_mtd_ttd$df==0,"z value","t value"), pname_mppo_mtd_ttd)
# mtests_mppo_mtd_ttd_exp_I = as.data.frame(mtests_mppo_mtd_ttd_exp_I)
# mtests_mppo_mtd_ttd_exp_I = mtests_mppo_mtd_ttd_exp_I[
#   which(mtests_mppo_mtd_ttd_exp_I$`Pr(>|z|)`<0.05),]

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMEANS
mixmodel.emm = emmeans(mixmodel, ~ MFPO * MTD * TTD)
contrast(mixmodel.emm, "eff", by = "MFPO")
contrast(mixmodel.emm, "eff", by = "MTD")
contrast(mixmodel.emm, "eff", by = "TTD")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans(mixmodel, c("MFPO","MTD","TTD"))

file_name = "interaction_plot_mfpo_mtd_ttd_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel, MTD ~ TTD | MFPO)
dev.off()

emm = emmeans(mixmodel, specs=pairwise~MFPO|MTD*TTD, type="response")

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_mfpo_mtd_ttd_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

# emmn = emmeans(mixmodel, ~MFPO|MTD*TTD, type="response")
# 
# file_name = "pwpp_mfpo_mtd_ttd_exp_I.pdf"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# pdf(file = figure_to_save, width=14, height=6)
# par(mar=c(4,5,0,0)+0.2)
# pwpp(emmn)
# dev.off()

# TABLES LATEX FOR SURVIE ~ MFPO * MTD * TTD
xtable(anova(mixmodel))
xtable(mtests_mppo_exp_I)
xtable(mtests_mtd_exp_I)
xtable(mtests_ttd_exp_I)
xtable(mtests_mppo_mtd_exp_I)
xtable(mtests_mppo_ttd_exp_I)
xtable(mtests_mtd_ttd_exp_I)
xtable(mtests_mppo_mtd_ttd_exp_I)
xtable(emm_df)

################################################################################
# MODELE DPRIME ~ UNCERTAINTY * SIMILARITY
################################################################################

setwd(exp_path)

dprime_Exp_I = read.table("psychoac_all_factors_Exp_I.csv", header=TRUE, 
                          sep=',', dec=".", fileEncoding="utf-8")

dprime_Exp_I = dprime_Exp_I[which(dprime_Exp_I$Sujet != 2 & 
                                    dprime_Exp_I$Sujet != 6 & 
                                    dprime_Exp_I$Sujet != 9),]

cond_uncertainty = c("29","115","463")
dprime_Exp_I$Uncertainty = factor(dprime_Exp_I$Uncertainty, 
                                  levels=cond_uncertainty)

cond_similarity = c("-80","-40","0","40","80")
dprime_Exp_I$Similarity = factor(dprime_Exp_I$Similarity, 
                                 levels=cond_similarity)

mixmodel = lme(dprime ~ Uncertainty * Similarity, data=dprime_Exp_I, 
               random=~1|Sujet, method="ML")

summary(mixmodel)
anova(mixmodel)

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# Uncertainty
mcp_uncertainty = glht(mixmodel, linfct=mcp(
  Uncertainty="Tukey", interaction_average=TRUE, covariate_average=TRUE))
summary(mcp_uncertainty, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_uncertainty)$test
mtests_uncertainty_exp_I=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_uncertainty$alternativ,
  less=paste("Pr(<", ifelse(mcp_uncertainty$df ==0, "z", "t"), ")", sep=""),
  greater=paste("Pr(>", ifelse(mcp_uncertainty$df==0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_uncertainty$df==0, "z", "t"),"|)",sep=""))
colnames(mtests_uncertainty_exp_I) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_uncertainty$df ==0, "z value", "t value"), pname)

# Similarity
mcp_similarity = glht(mixmodel, linfct=mcp(
  Similarity="Tukey", interaction_average=TRUE, covariate_average=TRUE))
summary(mcp_similarity, test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_similarity)$test
mtests_similarity_exp_I = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_similarity$alternativ,
  less = paste("Pr(<", ifelse(mcp_similarity$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_similarity$df == 0, "z", "t"), ")",sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_similarity$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_similarity_exp_I) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_similarity$df==0,"z value","t value"), pname)

# Uncertainty * Similarity
dprime_Exp_I$res = interaction(
  dprime_Exp_I$Uncertainty, dprime_Exp_I$Similarity)
dprime_Exp_I$res = factor(dprime_Exp_I$res)
model_uncertainty_similarity = lme(dprime ~ res, data=dprime_Exp_I, 
                                   random=~1|Sujet, method="ML")
mcp_unc_sim = glht(model_uncertainty_similarity, linfct=mcp(res="Tukey"))
summary(mcp_unc_sim, test=adjusted(type="bonferroni"))
pq_unc_sim = summary(mcp_unc_sim, test=adjusted(type="bonferroni"))$test
mtests_unc_sim_exp_I = cbind(pq_unc_sim$coefficients, pq_unc_sim$sigma, 
                              pq_unc_sim$tstat, pq_unc_sim$pvalues)
error_unc_sim = attr(pq_unc_sim$pvalues, "error")
pname_unc_sim = switch(
  mcp_unc_sim$alternativ,
  less = paste("Pr(<", ifelse(mcp_unc_sim$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_unc_sim$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_unc_sim$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_unc_sim_exp_I) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_unc_sim$df==0,"z value","t value"), pname_unc_sim)
# mtests_unc_sim_exp_I = as.data.frame(mtests_unc_sim_exp_I)
# mtests_unc_sim_exp_I = mtests_unc_sim_exp_I[
#   which(mtests_unc_sim_exp_I$`Pr(>|z|)`<0.05),]

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMEANS
mixmodel.emm = emmeans(mixmodel, ~ Uncertainty * Similarity)
contrast(mixmodel.emm, "eff", by = "Uncertainty")
contrast(mixmodel.emm, "eff", by = "Similarity")
contrast(mixmodel.emm, interaction = c("eff", "pairwise"), 
         adjust=multcomp_adjusted_type)

mixmodel.emm.interaction = emmeans(mixmodel, c("Uncertainty","Similarity"))

file_name = "interaction_plot_uncertainty_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel,  Uncertainty ~ Similarity)
dev.off()

file_name = "interaction_plot_similarity_uncertainty_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(mixmodel,  Similarity ~ Uncertainty)
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

# emmn = emmeans(mixmodel, ~Uncertainty|Similarity, type="response")
# 
# file_name = "pwpp_uncertainty_similarity_exp_I.pdf"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# pdf(file = figure_to_save, width=14, height=6)
# par(mar=c(4,5,0,0)+0.2)
# pwpp(emmn)
# dev.off()

# TABLES LATEX FOR SURVIE ~ Uncertainty * Similarity
xtable(anova(mixmodel))
xtable(mtests_uncertainty_exp_I)
xtable(mtests_similarity_exp_I)
xtable(mtests_unc_sim_exp_I)
xtable(emm_df)
