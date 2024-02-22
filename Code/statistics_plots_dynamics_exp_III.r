################################################################################
# SCRIPT FOR ANALYZING SURVIVAL DATA FROM BEHAVIORAL EXPERIMENT DURING PHD 
# THESIS FOR EXPERIMENT III
################################################################################

# CTRL + ALT + R

################################################################################
# LIBRARY LOAD
################################################################################

library(Matrix)
library(MASS)
library(bdsmatrix)
library(lattice)
library(mvtnorm)
library(TH.data)
library(cowplot)
library(survival)
library(coxme)
library(frailtyEM)
library(stringr)
library(multcomp)
library(xtable)
library(emmeans)
library(ggeffects)

setEPS()

################################################################################
# REFERENCING FOR PLOTTING
################################################################################

root_path = "/home/link/Documents/thèse_onera/"
article_path = "/home/link/Documents/thèse_onera/articles_alex/"
article_target = "article_1_PA/PloSOne/maxi_dossier/Exp_III/dynamics_stats"
slides_path = "/home/link/Documents/thèse_onera/diapos_phd_thesis/"
exp_path = "/home/link/Documents/thèse_onera/experimentation/psychophysique/"
slides_path_acc = "images/EEG/info_content"
slides_path_final = str_c(slides_path, slides_path_acc, sep="")

multcomp_adjusted_type = "fdr"
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")

################################################################################
# DATA LOAD
################################################################################

setwd("/home/link/Documents/thèse_onera/experimentation/psychophysique/")
setwd("data/test/testC/sujet_marseille/")
surv_df=read.table("surv.csv",header=TRUE,sep=',',dec=".",fileEncoding="utf-8")

################################################################################
# DATA CLEANING/RESHAPING
################################################################################

surv_df$RT = surv_df$RT/1000
surv_df$m_ppo = factor(surv_df$m_ppo)
surv_df$m_iti = factor(surv_df$m_iti)
surv_df$t_ra = factor(surv_df$t_ra)
surv_df$m_density = factor(surv_df$m_density)
surv_df$Uncertainty = factor(surv_df$Uncertainty)
surv_df = surv_df[which(surv_df$Sujet != 8),]
surv_df = surv_df[which(surv_df$Bloc != 1),] 
survie = Surv(surv_df$RT, surv_df$Surv)

# coln_original = c(
#   "X","Sujet","Bloc","Stim","Hits","FA","Miss","RC","RT","Diss","m_ppo","m_td",
#   "m_iti","t_pi","t_ra","t_td","m_density","Uncertainty","Surv")
# coln_modifed = c(
#   "X","Sujet","Bloc","Stim","Hits","FA","Miss","RC","RT","Similarity","FPO",
#   "MTD","MITI","TPI","TRA","TTD","Density","Uncertainty","Surv")
# 
# colnames(surv_df) = coln_modifed

################################################################################
# MODELE SURVIE ~ MPPO * MITI * TRA
################################################################################

# SURVIVAL DATA ANALYZING
model_nul = coxph(survie ~ m_ppo * m_iti * t_ra, data=surv_df)
model_mppo_miti_tra = coxph(survie ~ m_ppo * m_iti * t_ra + 
                             frailty(Sujet, distribution="gamma"), data=surv_df)
anova(model_nul, model_mppo_miti_tra)
summary(model_mppo_miti_tra)
anova(model_mppo_miti_tra)

# MODELS DIAGNOSTICS
mres = resid(model_mppo_miti_tra, type = "martingale")
csres = surv_df$Surv - mres
exp_dist = rexp(1000, rate=1)
par(mfrow=c(2,2))
hist(csres, breaks = 20, col=c("lightblue"), prob=TRUE, main="", 
     xlab="Cox-Snell Residuals", ylab="Density", ylim=c(0,1), xlim=c(0,6))
lines(density(csres), col="red")
hist(exp_dist, breaks = 20, col=c("lightblue"), prob=TRUE, main="", 
     xlab="Exponential Distribution", ylab="Density", ylim=c(0,1), xlim=c(0,6))
lines(density(exp_dist), col="red")
rsurv = survfit(Surv(csres, surv_df$Surv) ~ 1, type="fleming-harrington")
plot(0, 0, lty=1, type='n', xlim=c(0,3), ylim=c(0,3),
     main="Diagnostic plot of the frailty model in Experiment I",
     xlab="Cox-Snell residuals",
     ylab="Estimated Cumulative Hazard")
lines(rsurv$time, -log(rsurv$surv), type="s")
lines(c(0,3),c(0,3), col="red")

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# MPPO
mcp_mppo = glht(model_mppo_miti_tra, 
                linfct=mcp(m_ppo="Tukey", 
                           interaction_average = TRUE,
                           covariate_average = TRUE))
summary(mcp_mppo, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_mppo)$test
mtests_mppo_exp_III=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_mppo$alternativ,
  less=paste("Pr(<", ifelse(mcp_mppo$df ==0, "z", "t"), ")", sep = ""),
  greater=paste("Pr(>", ifelse(mcp_mppo$df==0, "z", "t"), ")", sep = ""),
  two.sided=paste("Pr(>|", ifelse(mcp_mppo$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_mppo_exp_III) = c("Estimate", "Std. Error",
                                 ifelse(mcp_mppo$df ==0, "z value", "t value"),
                                 pname)

# MITI
mcp_miti = glht(model_mppo_miti_tra, 
               linfct=mcp(m_iti="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_miti, test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_miti)$test
mtests_miti_exp_III = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_miti$alternativ,
  less = paste("Pr(<", ifelse(mcp_miti$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_miti$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_miti$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_miti_exp_III) = c("Estimate", "Std. Error",
                                ifelse(mcp_miti$df==0, "z value", "t value"),
                                pname)

# TRA
mcp_tra = glht(model_mppo_miti_tra, 
               linfct=mcp(t_ra="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_tra,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_tra)$test
mtests_tra_exp_III = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_tra$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_tra_exp_III) = c("Estimate", "Std. Error",
                                ifelse(mcp_tra$df==0, "z value", "t value"),
                                pname)

# MPPO * MITI
surv_df$res = interaction(surv_df$m_ppo, surv_df$m_iti)
surv_df$res = factor(surv_df$res)
model_mppo_miti = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_miti = glht(model_mppo_miti, linfct=mcp(res="Tukey"))
summary(mcp_mppo_miti, test=adjusted(type="bonferroni"))
pq_mppo_miti = summary(mcp_mppo_miti, test=adjusted(type="bonferroni"))$test
mtests_mppo_miti_exp_III = cbind(pq_mppo_miti$coefficients, pq_mppo_miti$sigma, 
                               pq_mppo_miti$tstat, pq_mppo_miti$pvalues)
error_mppo_miti = attr(pq_mppo_miti$pvalues, "error")
pname_mppo_miti = switch(
  mcp_mppo_miti$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_miti$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_miti$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_miti$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_miti_exp_III) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_miti$df==0,"z value","t value"), pname_mppo_miti)
mtests_mppo_miti_exp_III = as.data.frame(mtests_mppo_miti_exp_III)
mtests_mppo_miti_exp_III = mtests_mppo_miti_exp_III[
  which(mtests_mppo_miti_exp_III$`Pr(>|z|)`<0.05),]

# MPPO * TRA
surv_df$res = interaction(surv_df$m_ppo, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_mppo_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_tra = glht(model_mppo_tra, linfct=mcp(res="Tukey"))
summary(mcp_mppo_tra, test=adjusted(type="bonferroni"))
pq_mppo_tra = summary(mcp_mppo_tra, test=adjusted(type="bonferroni"))$test
mtests_mppo_tra_exp_III = cbind(pq_mppo_tra$coefficients, pq_mppo_tra$sigma, 
                               pq_mppo_tra$tstat, pq_mppo_tra$pvalues)
error_mppo_tra = attr(pq_mppo_tra$pvalues, "error")
pname_mppo_tra = switch(
  mcp_mppo_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_tra_exp_III) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_tra$df==0,"z value","t value"), pname_mppo_tra)
mtests_mppo_tra_exp_III = as.data.frame(mtests_mppo_tra_exp_III)
mtests_mppo_tra_exp_III = mtests_mppo_tra_exp_III[
  which(mtests_mppo_tra_exp_III$`Pr(>|z|)`<0.05),]

# MITI * TRA
surv_df$res = interaction(surv_df$m_iti, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_miti_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_miti_tra = glht(model_miti_tra, linfct=mcp(res="Tukey"))
summary(mcp_miti_tra, test=adjusted(type="bonferroni"))
pq_miti_tra = summary(mcp_miti_tra, test=adjusted(type="bonferroni"))$test
mtests_miti_tra_exp_III = cbind(pq_miti_tra$coefficients, pq_miti_tra$sigma, 
                              pq_miti_tra$tstat, pq_miti_tra$pvalues)
error_miti_tra = attr(pq_miti_tra$pvalues, "error")
pname_miti_tra = switch(
  mcp_miti_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_miti_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_miti_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_miti_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_miti_tra_exp_III) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_miti_tra$df==0,"z value","t value"), pname_miti_tra)
mtests_miti_tra_exp_III = as.data.frame(mtests_miti_tra_exp_III)
mtests_miti_tra_exp_III = mtests_miti_tra_exp_III[
  which(mtests_miti_tra_exp_III$`Pr(>|z|)`<0.05),]

# MPPO * MITI * TRA
surv_df$res = interaction(surv_df$m_ppo, surv_df$m_iti, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_mppo_miti_tra_res = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_miti_tra = glht(model_mppo_miti_tra_res, linfct=mcp(res="Tukey"))
summary(mcp_mppo_miti_tra, test=adjusted(type="bonferroni"))
pq_mppo_miti_tra = summary(mcp_mppo_miti_tra, 
                          test=adjusted(type="bonferroni"))$test
mtests_mppo_miti_tra_exp_III = cbind(
  pq_mppo_miti_tra$coefficients, pq_mppo_miti_tra$sigma, 
  pq_mppo_miti_tra$tstat, pq_mppo_miti_tra$pvalues)
error_mppo_miti_tra = attr(pq_mppo_miti_tra$pvalues, "error")
pname_mppo_miti_tra = switch(
  mcp_mppo_miti_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_miti_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_mppo_miti_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_mppo_miti_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_miti_tra_exp_III) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_miti_tra$df==0,"z value","t value"), pname_mppo_miti_tra)
mtests_mppo_miti_tra_exp_III = as.data.frame(mtests_mppo_miti_tra_exp_III)
mtests_mppo_miti_tra_exp_III = mtests_mppo_miti_tra_exp_III[
  which(mtests_mppo_miti_tra_exp_III$`Pr(>|z|)`<0.05),]

# # PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS
file_name = "interaction_plot_miti_tra_mfpo_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(model_mppo_miti_tra, m_iti ~ t_ra | m_ppo, CIs = TRUE)
dev.off()

emm = emmeans(model_mppo_miti_tra, specs=pairwise~t_ra|m_ppo*m_iti, 
              type="response", df=26)

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

# TABLES LATEX FOR SURVIE ~ MPPO * MITI * TRA
xtable(anova(model_nul, model_mppo_miti_tra))
xtable(anova(model_mppo_miti_tra))
xtable(mtests_mppo_exp_III)
xtable(mtests_miti_exp_III)
xtable(mtests_tra_exp_III)
xtable(mtests_mppo_miti_exp_III)
xtable(mtests_mppo_tra_exp_III)
xtable(mtests_miti_tra_exp_III)
xtable(mtests_mppo_miti_tra_exp_III)
xtable(emm_df)

# model_mppo_miti_tra = emmeans(model_mppo_miti_tra, ~ m_ppo * m_iti * t_ra)
# contrast(model_mppo_miti_tra, "eff", by = "m_ppo")
# contrast(model_mppo_miti_tra, "eff", by = "m_iti")
# contrast(model_mppo_miti_tra, "eff", by = "t_ra")
# contrast(model_mppo_miti_tra, interaction = c("eff", "pairwise"),
#          adjust = multcomp_adjusted_type)
# xtable(contrast(model_mppo_miti_tra, interaction = c("eff", "pairwise"),
#                 adjust = multcomp_adjusted_type))

################################################################################
# MODELE SURVIE ~ UNCERTAINTY * TRA
################################################################################

# SURVIVAL DATA ANALYZING
model_nul = coxph(survie ~ Uncertainty * t_ra, data=surv_df)
model_unc_tra = coxph(survie ~ Uncertainty * t_ra + 
                             frailty(Sujet, distribution="gamma"), data=surv_df)
anova(model_nul, model_unc_tra)
summary(model_unc_tra)
anova(model_unc_tra)

# MODELS DIAGNOSTICS
mres = resid(model_unc_tra, type = "martingale")
csres = surv_df$Surv - mres
exp_dist = rexp(1000, rate=1)
par(mfrow=c(2,2))
hist(csres, breaks = 20, col=c("lightblue"), prob=TRUE, main="", 
     xlab="Cox-Snell Residuals", ylab="Density", ylim=c(0,1), xlim=c(0,6))
lines(density(csres), col="red")
hist(exp_dist, breaks = 20, col=c("lightblue"), prob=TRUE, main="", 
     xlab="Exponential Distribution", ylab="Density", ylim=c(0,1), xlim=c(0,6))
lines(density(exp_dist), col="red")
rsurv = survfit(Surv(csres, surv_df$Surv) ~ 1, type="fleming-harrington")
plot(0, 0, lty=1, type='n', xlim=c(0,3), ylim=c(0,3),
     main="Diagnostic plot of the frailty model in Experiment I",
     xlab="Cox-Snell residuals",
     ylab="Estimated Cumulative Hazard")
lines(rsurv$time, -log(rsurv$surv), type="s")
lines(c(0,3),c(0,3), col="red")

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# UNCERTAINTY
mcp_unc = glht(model_unc_tra, 
               linfct=mcp(Uncertainty="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_unc, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_unc)$test
mtests_unc_exp_III=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_unc$alternativ,
  less=paste("Pr(<", ifelse(mcp_unc$df ==0, "z", "t"), ")", sep = ""),
  greater=paste("Pr(>", ifelse(mcp_unc$df==0, "z", "t"), ")", sep = ""),
  two.sided=paste("Pr(>|", ifelse(mcp_unc$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_unc_exp_III) = c("Estimate", "Std. Error",
                                ifelse(mcp_unc$df ==0, "z value", "t value"),
                                pname)

# TRA
mcp_tra = glht(model_unc_tra, 
               linfct=mcp(t_ra="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_tra,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_tra)$test
mtests_tra_exp_III = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_tra$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_tra_exp_III) = c("Estimate", "Std. Error",
                                ifelse(mcp_tra$df==0, "z value", "t value"),
                                pname)

# UNCERTAINTY * TRA
surv_df$res = interaction(surv_df$Uncertainty, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_unc_tra_res = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_unc_tra = glht(model_unc_tra_res, linfct=mcp(res="Tukey"))
summary(mcp_unc_tra, test=adjusted(type="bonferroni"))
pq_unc_tra = summary(mcp_unc_tra, test=adjusted(type="bonferroni"))$test
mtests_unc_tra_exp_III = cbind(pq_unc_tra$coefficients, pq_unc_tra$sigma, 
                              pq_unc_tra$tstat, pq_unc_tra$pvalues)
error_unc_tra = attr(pq_unc_tra$pvalues, "error")
pname_unc_tra = switch(
  mcp_unc_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_unc_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_unc_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_unc_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_unc_tra_exp_III) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_unc_tra$df==0,"z value","t value"), pname_unc_tra)
mtests_unc_tra_exp_III = as.data.frame(mtests_unc_tra_exp_III)
mtests_unc_tra_exp_III = mtests_unc_tra_exp_III[
  which(mtests_unc_tra_exp_III$`Pr(>|z|)`<0.05),]

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS
file_name = "interaction_plot_targetrate_uncertainty_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(model_unc_tra, t_ra ~ Uncertainty)
dev.off()

emm = emmeans(model_unc_tra, specs=pairwise~t_ra|Uncertainty, type="response", 
              df=26)

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_targetrate_uncertainty_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emmeans(model_unc_tra, specs=pairwise~Uncertainty|t_ra, 
             type="response", df=26), comparisons = TRUE)
dev.off()

# TABLES LATEX FOR SURVIE ~ UNCERTAINTY * SIMILARITY * TRA
xtable(anova(model_nul, model_unc_tra))
xtable(anova(model_unc_tra))
xtable(mtests_unc_exp_III)
xtable(mtests_tra_exp_III)
xtable(mtests_unc_tra_exp_III)
xtable(emm_df)

# model_unc_tra = emmeans(model_unc_tra, ~ Uncertainty * t_ra)
# contrast(model_unc_tra, "eff", by = "Uncertainty")
# contrast(model_unc_tra, "eff", by = "t_ra")
# contrast(model_unc_tra, interaction = c("eff", "pairwise"),
#          adjust = multcomp_adjusted_type)
# xtable(contrast(model_unc_tra, interaction = c("eff", "pairwise"),
#                 adjust = multcomp_adjusted_type))

################################################################################
# MODELE SURVIE ~ DENSITY * TRA
################################################################################

# SURVIVAL DATA ANALYZING
model_nul = coxph(survie ~ m_density * t_ra, data=surv_df)
model_den_tra = coxph(survie ~ m_density * t_ra + 
                        frailty(Sujet, distribution="gamma"), data=surv_df)
anova(model_nul, model_den_tra)
summary(model_den_tra)
anova(model_den_tra)

# MODELS DIAGNOSTICS
mres = resid(model_den_tra, type = "martingale")
csres = surv_df$Surv - mres
exp_dist = rexp(1000, rate=1)
par(mfrow=c(2,2))
hist(csres, breaks = 20, col=c("lightblue"), prob=TRUE, main="", 
     xlab="Cox-Snell Residuals", ylab="Density", ylim=c(0,1), xlim=c(0,6))
lines(density(csres), col="red")
hist(exp_dist, breaks = 20, col=c("lightblue"), prob=TRUE, main="", 
     xlab="Exponential Distribution", ylab="Density", ylim=c(0,1), xlim=c(0,6))
lines(density(exp_dist), col="red")
rsurv = survfit(Surv(csres, surv_df$Surv) ~ 1, type="fleming-harrington")
plot(0, 0, lty=1, type='n', xlim=c(0,3), ylim=c(0,3),
     main="Diagnostic plot of the frailty model in Experiment I",
     xlab="Cox-Snell residuals",
     ylab="Estimated Cumulative Hazard")
lines(rsurv$time, -log(rsurv$surv), type="s")
lines(c(0,3),c(0,3), col="red")

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# DENSITY
mcp_den = glht(model_den_tra, 
               linfct=mcp(m_density="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_den, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_den)$test
mtests_den_exp_III=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_den$alternativ,
  less=paste("Pr(<", ifelse(mcp_den$df ==0, "z", "t"), ")", sep = ""),
  greater=paste("Pr(>", ifelse(mcp_den$df==0, "z", "t"), ")", sep = ""),
  two.sided=paste("Pr(>|", ifelse(mcp_den$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_den_exp_III) = c("Estimate", "Std. Error",
                                ifelse(mcp_den$df ==0, "z value", "t value"),
                                pname)

# TRA
mcp_tra = glht(model_den_tra, 
               linfct=mcp(t_ra="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_tra,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_tra)$test
mtests_tra_exp_III = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_tra$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_tra_exp_III) = c("Estimate", "Std. Error",
                                ifelse(mcp_tra$df==0, "z value", "t value"),
                                pname)

# DENSITY * TRA
surv_df$res = interaction(surv_df$m_density, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_den_tra_res = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_den_tra = glht(model_den_tra_res, linfct=mcp(res="Tukey"))
summary(mcp_den_tra, test=adjusted(type="bonferroni"))
pq_den_tra = summary(mcp_den_tra, test=adjusted(type="bonferroni"))$test
mtests_den_tra_exp_III = cbind(pq_den_tra$coefficients, pq_den_tra$sigma, 
                              pq_den_tra$tstat, pq_den_tra$pvalues)
error_den_tra = attr(pq_den_tra$pvalues, "error")
pname_den_tra = switch(
  mcp_den_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_den_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_den_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_den_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_den_tra_exp_III) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_den_tra$df==0,"z value","t value"), pname_den_tra)
mtests_den_tra_exp_III = as.data.frame(mtests_den_tra_exp_III)
mtests_den_tra_exp_III = mtests_den_tra_exp_III[
  which(mtests_den_tra_exp_III$`Pr(>|z|)`<0.05),]

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS
file_name = "interaction_plot_density_targetrate_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(model_den_tra, t_ra ~ m_density)
dev.off()

emm = emmeans(model_den_tra, specs=pairwise~t_ra|m_density, type="response", 
              df=26)

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_density_targetrate_exp_III.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emmeans(model_den_tra, specs=pairwise~m_density|t_ra, type="response", 
             df=26), comparisons = TRUE)
dev.off()

# TABLES LATEX FOR SURVIE ~ DENSITY * SIMILARITY * TRA
xtable(anova(model_nul, model_den_tra))
xtable(anova(model_den_tra))
xtable(mtests_den_exp_III)
xtable(mtests_tra_exp_III)
xtable(mtests_den_tra_exp_III)
xtable(emm_df)

# model_den_tra = emmeans(model_den_tra, ~ m_density * t_ra)
# contrast(model_den_tra, "eff", by = "m_density")
# contrast(model_den_tra, "eff", by = "t_ra")
# contrast(model_den_tra, interaction = c("eff", "pairwise"),
#          adjust = multcomp_adjusted_type)
# xtable(contrast(model_den_tra, interaction = c("eff", "pairwise"),
#                 adjust = multcomp_adjusted_type))

################################################################################
# CUMULATIVE DISTRIBUTION FUNCTION EN COULEUR ET TRAIT POINTILLE
################################################################################

# # survie ~ m_ppo
# c_mppo = survfit(survie ~ m_ppo, data = surv_df)
# file_name = "cdf_mppo_exp_III.eps"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# postscript(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(c_mppo, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3), 
#      main = "", xlab = "Time (sec)", ylab = "F(t)", 
#      las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
# legend("bottomright", c("4fpo", "32fpo", "64fpo"), lty=c(1,2,3), lwd=c(3,3,3), 
#        col = c("green", "red", "blue"), cex=1.6)
# dev.off()
# 
# # survie ~ m_iti
# c_miti = survfit(survie ~ m_iti, data = surv_df)
# file_name = "cdf_miti_exp_III.eps"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# postscript(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(c_miti, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
#      main = "", xlab = "Time (sec)", ylab = "F(t)", 
#      las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
# legend("bottomright", c("200ms","600ms","1200ms"), lty=c(1,2,3), 
#        lwd=c(3,3,3), col = c("green", "red", "blue"), cex=1.6)
# dev.off()
# 
# # survie ~ t_ra
# c_tra = survfit(survie ~ t_ra, data = surv_df)
# file_name = "cdf_tra_exp_III.eps"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# postscript(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(c_tra, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
#      main = "", xlab = "Time (sec)", ylab = "F(t)", 
#      las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
# legend("bottomright", c("1Hz","2Hz","5Hz"), lty=c(1,2,3), 
#        lwd=c(3,3,3), col = c("green", "red", "blue"), cex=1.6)
# dev.off()
# 
# # survie ~ Uncertainty
# c_uncertainty = survfit(survie ~ Uncertainty, data = surv_df)
# file_name = "cdf_uncertainty_exp_III.eps"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# postscript(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(c_uncertainty, 
#      col = c("chartreuse4","green2","darkgreen","red1","red3","darkred",
#              "deepskyblue1","blue1","blue4"), 
#      lwd=c(3,3,3), main = "", xlab = "Time (sec)", ylab = "F(t)", las=1, 
#      cex.axis=1.6, cex.lab=1.6, fun="F")
# legend("bottomright", 
#        col = c("chartreuse4","green2","darkgreen","red1","red3","darkred",
#                "deepskyblue1","blue1","blue4"), 
#        c("84","110","123","169", "221", "246", "339", "442", "492"),
#        lwd=c(3,3,3), cex=1.2)
# dev.off()
# 
