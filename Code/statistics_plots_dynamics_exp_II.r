################################################################################
# SCRIPT FOR ANALYZING SURVIVAL DATA FROM BEHAVIORAL EXPERIMENT DURING PHD 
# THESIS FOR EXPERIMENT II
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
article_target = "article_1_PA/PloSOne/maxi_dossier/Exp_II/dynamics_stats"
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
setwd("data/test/testB/sujet_marseille/")
surv_df=read.table("surv.csv",header=TRUE,sep=',',dec=".",fileEncoding="utf-8")

################################################################################
# DATA CLEANING/RESHAPING
################################################################################

surv_df$RT = surv_df$RT/1000
surv_df$m_ppo = factor(surv_df$m_ppo)
surv_df$m_td = factor(surv_df$m_td)
surv_df$t_ra = factor(surv_df$t_ra)
surv_df$Diss = factor(surv_df$Diss)
surv_df$m_density = factor(surv_df$m_density)
surv_df$Uncertainty = factor(surv_df$Uncertainty)
surv_df = surv_df[which(surv_df$Sujet != 1 & surv_df$Sujet != 2 
                        & surv_df$Sujet != 6 & surv_df$Sujet != 8),]
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
# MODELE SURVIE ~ MPPO * MTD * TRA
################################################################################

# SURVIVAL DATA ANALYZING
model_nul = coxph(survie ~ m_ppo * m_td * t_ra, data=surv_df)
model_mppo_mtd_tra = coxph(survie ~ m_ppo * m_td * t_ra + 
                             frailty(Sujet, distribution="gamma"), data=surv_df)
anova(model_nul, model_mppo_mtd_tra)
summary(model_mppo_mtd_tra)
anova(model_mppo_mtd_tra)

# MODELS DIAGNOSTICS
mres = resid(model_mppo_mtd_tra, type = "martingale")
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
mcp_mppo = glht(model_mppo_mtd_tra, 
                linfct=mcp(m_ppo="Tukey", 
                           interaction_average = TRUE,
                           covariate_average = TRUE))
summary(mcp_mppo, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_mppo)$test
mtests_mppo_exp_II=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_mppo$alternativ,
  less=paste("Pr(<", ifelse(mcp_mppo$df ==0, "z", "t"), ")", sep = ""),
  greater=paste("Pr(>", ifelse(mcp_mppo$df==0, "z", "t"), ")", sep = ""),
  two.sided=paste("Pr(>|", ifelse(mcp_mppo$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_mppo_exp_II) = c("Estimate", "Std. Error",
                                ifelse(mcp_mppo$df ==0, "z value", "t value"),
                                pname)

# MTD
mcp_mtd = glht(model_mppo_mtd_tra, 
               linfct=mcp(m_td="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_mtd, test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_mtd)$test
mtests_mtd_exp_II = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_mtd$alternativ,
  less = paste("Pr(<", ifelse(mcp_mtd$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_mtd$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_mtd$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_mtd_exp_II) = c("Estimate", "Std. Error",
                               ifelse(mcp_mtd$df==0, "z value", "t value"),
                               pname)

# TRA
mcp_tra = glht(model_mppo_mtd_tra, 
               linfct=mcp(t_ra="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_tra,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_tra)$test
mtests_tra_exp_II = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_tra$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_tra_exp_II) = c("Estimate", "Std. Error",
                               ifelse(mcp_tra$df==0, "z value", "t value"),
                               pname)

# MPPO * MTD
surv_df$res = interaction(surv_df$m_ppo, surv_df$m_td)
surv_df$res = factor(surv_df$res)
model_mppo_mtd = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_mtd = glht(model_mppo_mtd, linfct=mcp(res="Tukey"))
summary(mcp_mppo_mtd, test=adjusted(type="bonferroni"))
pq_mppo_mtd = summary(mcp_mppo_mtd, test=adjusted(type="bonferroni"))$test
mtests_mppo_mtd_exp_II = cbind(pq_mppo_mtd$coefficients, pq_mppo_mtd$sigma, 
                              pq_mppo_mtd$tstat, pq_mppo_mtd$pvalues)
error_mppo_mtd = attr(pq_mppo_mtd$pvalues, "error")
pname_mppo_mtd = switch(
  mcp_mppo_mtd$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_mtd$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_mtd$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_mtd$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_mtd_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_mtd$df==0,"z value","t value"), pname_mppo_mtd)
mtests_mppo_mtd_exp_II = as.data.frame(mtests_mppo_mtd_exp_II)
mtests_mppo_mtd_exp_II = mtests_mppo_mtd_exp_II[
  which(mtests_mppo_mtd_exp_II$`Pr(>|z|)`<0.05),]

# MPPO * TRA
surv_df$res = interaction(surv_df$m_ppo, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_mppo_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_tra = glht(model_mppo_tra, linfct=mcp(res="Tukey"))
summary(mcp_mppo_tra, test=adjusted(type="bonferroni"))
pq_mppo_tra = summary(mcp_mppo_tra, test=adjusted(type="bonferroni"))$test
mtests_mppo_tra_exp_II = cbind(pq_mppo_tra$coefficients, pq_mppo_tra$sigma, 
                              pq_mppo_tra$tstat, pq_mppo_tra$pvalues)
error_mppo_tra = attr(pq_mppo_tra$pvalues, "error")
pname_mppo_tra = switch(
  mcp_mppo_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_tra$df==0,"z value","t value"), pname_mppo_tra)
mtests_mppo_tra_exp_II = as.data.frame(mtests_mppo_tra_exp_II)
mtests_mppo_tra_exp_II = mtests_mppo_tra_exp_II[
  which(mtests_mppo_tra_exp_II$`Pr(>|z|)`<0.05),]

# MTD * TRA
surv_df$res = interaction(surv_df$m_td, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_mtd_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mtd_tra = glht(model_mtd_tra, linfct=mcp(res="Tukey"))
summary(mcp_mtd_tra, test=adjusted(type="bonferroni"))
pq_mtd_tra = summary(mcp_mtd_tra, test=adjusted(type="bonferroni"))$test
mtests_mtd_tra_exp_II = cbind(pq_mtd_tra$coefficients, pq_mtd_tra$sigma, 
                             pq_mtd_tra$tstat, pq_mtd_tra$pvalues)
error_mtd_tra = attr(pq_mtd_tra$pvalues, "error")
pname_mtd_tra = switch(
  mcp_mtd_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_mtd_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mtd_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mtd_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mtd_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mtd_tra$df==0,"z value","t value"), pname_mtd_tra)
mtests_mtd_tra_exp_II = as.data.frame(mtests_mtd_tra_exp_II)
mtests_mtd_tra_exp_II = mtests_mtd_tra_exp_II[
  which(mtests_mtd_tra_exp_II$`Pr(>|z|)`<0.05),]

# MPPO * MTD * TRA
surv_df$res = interaction(surv_df$m_ppo, surv_df$m_td, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_mppo_mtd_tra_res = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_mtd_tra = glht(model_mppo_mtd_tra_res, linfct=mcp(res="Tukey"))
summary(mcp_mppo_mtd_tra, test=adjusted(type="bonferroni"))
pq_mppo_mtd_tra = summary(mcp_mppo_mtd_tra, 
                          test=adjusted(type="bonferroni"))$test
mtests_mppo_mtd_tra_exp_II = cbind(
  pq_mppo_mtd_tra$coefficients, pq_mppo_mtd_tra$sigma, 
  pq_mppo_mtd_tra$tstat, pq_mppo_mtd_tra$pvalues)
error_mppo_mtd_tra = attr(pq_mppo_mtd_tra$pvalues, "error")
pname_mppo_mtd_tra = switch(
  mcp_mppo_mtd_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_mtd_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_mppo_mtd_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_mppo_mtd_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_mtd_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_mtd_tra$df==0,"z value","t value"), pname_mppo_mtd_tra)
mtests_mppo_mtd_tra_exp_II = as.data.frame(mtests_mppo_mtd_tra_exp_II)
mtests_mppo_mtd_tra_exp_II = mtests_mppo_mtd_tra_exp_II[
  which(mtests_mppo_mtd_tra_exp_II$`Pr(>|z|)`<0.05),]

# # PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS

file_name = "interaction_plot_mtd_tra_mfpo_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(model_mppo_mtd_tra, m_td ~ t_ra | m_ppo, CIs = TRUE)
dev.off()

emm = emmeans(model_mppo_mtd_tra, specs=pairwise~t_ra|m_ppo*m_td, 
              type="response", df=26)

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

# TABLES LATEX FOR SURVIE ~ MPPO * MTD * TRA
xtable(anova(model_nul, model_mppo_mtd_tra))
xtable(anova(model_mppo_mtd_tra))
xtable(mtests_mppo_exp_II)
xtable(mtests_mtd_exp_II)
xtable(mtests_tra_exp_II)
xtable(mtests_mppo_mtd_exp_II)
xtable(mtests_mppo_tra_exp_II)
xtable(mtests_mtd_tra_exp_II)
xtable(mtests_mppo_mtd_tra_exp_II)
xtable(emm_df)

# model_mppo_mtd_tra = emmeans(model_mppo_mtd_tra, ~ m_ppo * m_td * t_ra)
# contrast(model_mppo_mtd_tra, "eff", by = "m_ppo")
# contrast(model_mppo_mtd_tra, "eff", by = "m_td")
# contrast(model_mppo_mtd_tra, "eff", by = "t_ra")
# contrast(model_mppo_mtd_tra, interaction = c("eff", "pairwise"),
#          adjust = multcomp_adjusted_type)
# xtable(contrast(model_mppo_mtd_tra, interaction = c("eff", "pairwise"),
#                 adjust = multcomp_adjusted_type))

################################################################################
# MODELE SURVIE ~ MPPO * SIMILARITY * TRA
################################################################################

# SURVIVAL DATA ANALYZING
model_nul = coxph(survie ~ m_ppo * Diss * t_ra, data=surv_df)
model_mppo_diss_tra = coxph(survie ~ m_ppo * Diss * t_ra + 
                             frailty(Sujet, distribution="gamma"), data=surv_df)
anova(model_nul, model_mppo_diss_tra)
summary(model_mppo_diss_tra)
anova(model_mppo_diss_tra)

# MODELS DIAGNOSTICS
mres = resid(model_mppo_diss_tra, type = "martingale")
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
mcp_mppo = glht(model_mppo_diss_tra, 
                linfct=mcp(m_ppo="Tukey", 
                           interaction_average = TRUE,
                           covariate_average = TRUE))
summary(mcp_mppo, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_mppo)$test
mtests_mppo_exp_II=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_mppo$alternativ,
  less=paste("Pr(<", ifelse(mcp_mppo$df ==0, "z", "t"), ")", sep = ""),
  greater=paste("Pr(>", ifelse(mcp_mppo$df==0, "z", "t"), ")", sep = ""),
  two.sided=paste("Pr(>|", ifelse(mcp_mppo$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_mppo_exp_II) = c("Estimate", "Std. Error",
                                 ifelse(mcp_mppo$df ==0, "z value", "t value"),
                                 pname)

# SIMILARITY
mcp_diss = glht(model_mppo_diss_tra, 
               linfct=mcp(Diss="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_diss, test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_diss)$test
mtests_diss_exp_II = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_diss$alternativ,
  less = paste("Pr(<", ifelse(mcp_diss$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_diss$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_diss$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_diss_exp_II) = c("Estimate", "Std. Error",
                                ifelse(mcp_diss$df==0, "z value", "t value"),
                                pname)

# TRA
mcp_tra = glht(model_mppo_diss_tra, 
               linfct=mcp(t_ra="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_tra,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_tra)$test
mtests_tra_exp_II = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_tra$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_tra_exp_II) = c("Estimate", "Std. Error",
                                ifelse(mcp_tra$df==0, "z value", "t value"),
                                pname)

# MPPO * SIMILARITY
surv_df$res = interaction(surv_df$m_ppo, surv_df$Diss)
surv_df$res = factor(surv_df$res)
model_mppo_diss = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_diss = glht(model_mppo_diss, linfct=mcp(res="Tukey"))
summary(mcp_mppo_diss, test=adjusted(type="bonferroni"))
pq_mppo_diss = summary(mcp_mppo_diss, test=adjusted(type="bonferroni"))$test
mtests_mppo_diss_exp_II = cbind(pq_mppo_diss$coefficients, pq_mppo_diss$sigma, 
                               pq_mppo_diss$tstat, pq_mppo_diss$pvalues)
error_mppo_diss = attr(pq_mppo_diss$pvalues, "error")
pname_mppo_diss = switch(
  mcp_mppo_diss$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_diss$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_diss$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_diss$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_diss_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_diss$df==0,"z value","t value"), pname_mppo_diss)
mtests_mppo_diss_exp_II = as.data.frame(mtests_mppo_diss_exp_II)
mtests_mppo_diss_exp_II = mtests_mppo_diss_exp_II[
  which(mtests_mppo_diss_exp_II$`Pr(>|z|)`<0.05),]

# MPPO * TRA
surv_df$res = interaction(surv_df$m_ppo, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_mppo_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_tra = glht(model_mppo_tra, linfct=mcp(res="Tukey"))
summary(mcp_mppo_tra, test=adjusted(type="bonferroni"))
pq_mppo_tra = summary(mcp_mppo_tra, test=adjusted(type="bonferroni"))$test
mtests_mppo_tra_exp_II = cbind(pq_mppo_tra$coefficients, pq_mppo_tra$sigma, 
                               pq_mppo_tra$tstat, pq_mppo_tra$pvalues)
error_mppo_tra = attr(pq_mppo_tra$pvalues, "error")
pname_mppo_tra = switch(
  mcp_mppo_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_tra$df==0,"z value","t value"), pname_mppo_tra)
mtests_mppo_tra_exp_II = as.data.frame(mtests_mppo_tra_exp_II)
mtests_mppo_tra_exp_II = mtests_mppo_tra_exp_II[
  which(mtests_mppo_tra_exp_II$`Pr(>|z|)`<0.05),]

# SIMILARITY * TRA
surv_df$res = interaction(surv_df$Diss, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_diss_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_diss_tra = glht(model_diss_tra, linfct=mcp(res="Tukey"))
summary(mcp_diss_tra, test=adjusted(type="bonferroni"))
pq_diss_tra = summary(mcp_diss_tra, test=adjusted(type="bonferroni"))$test
mtests_diss_tra_exp_II = cbind(pq_diss_tra$coefficients, pq_diss_tra$sigma, 
                              pq_diss_tra$tstat, pq_diss_tra$pvalues)
error_diss_tra = attr(pq_diss_tra$pvalues, "error")
pname_diss_tra = switch(
  mcp_diss_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_diss_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_diss_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_diss_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_diss_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_diss_tra$df==0,"z value","t value"), pname_diss_tra)
mtests_diss_tra_exp_II = as.data.frame(mtests_diss_tra_exp_II)
mtests_diss_tra_exp_II = mtests_diss_tra_exp_II[
  which(mtests_diss_tra_exp_II$`Pr(>|z|)`<0.05),]

# MPPO * SIMILARITY * TRA
surv_df$res = interaction(surv_df$m_ppo, surv_df$Diss, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_mppo_diss_tra_res = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_diss_tra = glht(model_mppo_diss_tra_res, linfct=mcp(res="Tukey"))
summary(mcp_mppo_diss_tra, test=adjusted(type="bonferroni"))
pq_mppo_diss_tra = summary(mcp_mppo_diss_tra, 
                          test=adjusted(type="bonferroni"))$test
mtests_mppo_diss_tra_exp_II = cbind(
  pq_mppo_diss_tra$coefficients, pq_mppo_diss_tra$sigma, 
  pq_mppo_diss_tra$tstat, pq_mppo_diss_tra$pvalues)
error_mppo_diss_tra = attr(pq_mppo_diss_tra$pvalues, "error")
pname_mppo_diss_tra = switch(
  mcp_mppo_diss_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_diss_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_mppo_diss_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_mppo_diss_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_diss_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_mppo_diss_tra$df==0,"z value","t value"), pname_mppo_diss_tra)
mtests_mppo_diss_tra_exp_II = as.data.frame(mtests_mppo_diss_tra_exp_II)
mtests_mppo_diss_tra_exp_II = mtests_mppo_diss_tra_exp_II[
  which(mtests_mppo_diss_tra_exp_II$`Pr(>|z|)`<0.05),]

# # PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS

file_name = "interaction_plot_sim_tra_mfpo_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(model_mppo_diss_tra, Diss ~ t_ra | m_ppo, CIs = TRUE)
dev.off()

emm = emmeans(model_mppo_diss_tra, specs=pairwise~t_ra|m_ppo*Diss, 
              type="response", df=26)

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

# TABLES LATEX FOR SURVIE ~ MPPO * SIMILARITY * TRA
xtable(anova(model_nul, model_mppo_diss_tra))
xtable(anova(model_mppo_diss_tra))
xtable(mtests_mppo_exp_II)
xtable(mtests_diss_exp_II)
xtable(mtests_tra_exp_II)
xtable(mtests_mppo_diss_exp_II)
xtable(mtests_mppo_tra_exp_II)
xtable(mtests_diss_tra_exp_II)
xtable(mtests_mppo_diss_tra_exp_II)
xtable(emm_df)

# model_mppo_diss_tra = emmeans(model_mppo_diss_tra, ~ m_ppo * Diss * t_ra)
# contrast(model_mppo_diss_tra, "eff", by = "m_ppo")
# contrast(model_mppo_diss_tra, "eff", by = "Diss")
# contrast(model_mppo_diss_tra, "eff", by = "t_ra")
# contrast(model_mppo_diss_tra, interaction = c("eff", "pairwise"),
#          adjust = multcomp_adjusted_type)
# xtable(contrast(model_mppo_diss_tra, interaction = c("eff", "pairwise"),
#                 adjust = multcomp_adjusted_type))

################################################################################
# MODELE SURVIE ~ UNCERTAINTY * SIMILARITY * TRA
################################################################################

# SURVIVAL DATA ANALYZING
model_nul = coxph(survie ~ Uncertainty * Diss * t_ra, data=surv_df)
model_unc_diss_tra = coxph(survie ~ Uncertainty * Diss * t_ra + 
                             frailty(Sujet, distribution="gamma"), data=surv_df)
anova(model_nul, model_unc_diss_tra)
summary(model_unc_diss_tra)
anova(model_unc_diss_tra)

# MODELS DIAGNOSTICS
mres = resid(model_unc_diss_tra, type = "martingale")
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
mcp_unc = glht(model_unc_diss_tra, 
                linfct=mcp(Uncertainty="Tukey", 
                           interaction_average = TRUE,
                           covariate_average = TRUE))
summary(mcp_unc, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_unc)$test
mtests_unc_exp_II=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_unc$alternativ,
  less=paste("Pr(<", ifelse(mcp_unc$df ==0, "z", "t"), ")", sep = ""),
  greater=paste("Pr(>", ifelse(mcp_unc$df==0, "z", "t"), ")", sep = ""),
  two.sided=paste("Pr(>|", ifelse(mcp_unc$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_unc_exp_II) = c("Estimate", "Std. Error",
                                 ifelse(mcp_unc$df ==0, "z value", "t value"),
                                 pname)

# SIMILARITY
mcp_diss = glht(model_unc_diss_tra, 
                linfct=mcp(Diss="Tukey", 
                           interaction_average = TRUE,
                           covariate_average = TRUE))
summary(mcp_diss, test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_diss)$test
mtests_diss_exp_II = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_diss$alternativ,
  less = paste("Pr(<", ifelse(mcp_diss$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_diss$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_diss$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_diss_exp_II) = c("Estimate", "Std. Error",
                                 ifelse(mcp_diss$df==0, "z value", "t value"),
                                 pname)

# TRA
mcp_tra = glht(model_unc_diss_tra, 
               linfct=mcp(t_ra="Tukey", 
                          interaction_average = TRUE,
                          covariate_average = TRUE))
summary(mcp_tra,test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_tra)$test
mtests_tra_exp_II = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_tra$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_tra_exp_II) = c("Estimate", "Std. Error",
                                ifelse(mcp_tra$df==0, "z value", "t value"),
                                pname)

# UNCERTAINTY * SIMILARITY
surv_df$res = interaction(surv_df$Uncertainty, surv_df$Diss)
surv_df$res = factor(surv_df$res)
model_unc_diss = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_unc_diss = glht(model_unc_diss, linfct=mcp(res="Tukey"))
summary(mcp_unc_diss, test=adjusted(type="bonferroni"))
pq_unc_diss = summary(mcp_unc_diss, test=adjusted(type="bonferroni"))$test
mtests_unc_diss_exp_II = cbind(pq_unc_diss$coefficients, pq_unc_diss$sigma, 
                                pq_unc_diss$tstat, pq_unc_diss$pvalues)
error_unc_diss = attr(pq_unc_diss$pvalues, "error")
pname_unc_diss = switch(
  mcp_unc_diss$alternativ,
  less = paste("Pr(<", ifelse(mcp_unc_diss$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_unc_diss$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_unc_diss$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_unc_diss_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_unc_diss$df==0,"z value","t value"), pname_unc_diss)
mtests_unc_diss_exp_II = as.data.frame(mtests_unc_diss_exp_II)
mtests_unc_diss_exp_II = mtests_unc_diss_exp_II[
  which(mtests_unc_diss_exp_II$`Pr(>|z|)`<0.05),]

# UNCERTAINTY * TRA
surv_df$res = interaction(surv_df$Uncertainty, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_unc_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_unc_tra = glht(model_unc_tra, linfct=mcp(res="Tukey"))
summary(mcp_unc_tra, test=adjusted(type="bonferroni"))
pq_unc_tra = summary(mcp_unc_tra, test=adjusted(type="bonferroni"))$test
mtests_unc_tra_exp_II = cbind(pq_unc_tra$coefficients, pq_unc_tra$sigma, 
                               pq_unc_tra$tstat, pq_unc_tra$pvalues)
error_unc_tra = attr(pq_unc_tra$pvalues, "error")
pname_unc_tra = switch(
  mcp_unc_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_unc_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_unc_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_unc_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_unc_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_unc_tra$df==0,"z value","t value"), pname_unc_tra)
mtests_unc_tra_exp_II = as.data.frame(mtests_unc_tra_exp_II)
mtests_unc_tra_exp_II = mtests_unc_tra_exp_II[
  which(mtests_unc_tra_exp_II$`Pr(>|z|)`<0.05),]

# SIMILARITY * TRA
surv_df$res = interaction(surv_df$Diss, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_diss_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_diss_tra = glht(model_diss_tra, linfct=mcp(res="Tukey"))
summary(mcp_diss_tra, test=adjusted(type="bonferroni"))
pq_diss_tra = summary(mcp_diss_tra, test=adjusted(type="bonferroni"))$test
mtests_diss_tra_exp_II = cbind(pq_diss_tra$coefficients, pq_diss_tra$sigma, 
                               pq_diss_tra$tstat, pq_diss_tra$pvalues)
error_diss_tra = attr(pq_diss_tra$pvalues, "error")
pname_diss_tra = switch(
  mcp_diss_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_diss_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_diss_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_diss_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_diss_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_diss_tra$df==0,"z value","t value"), pname_diss_tra)
mtests_diss_tra_exp_II = as.data.frame(mtests_diss_tra_exp_II)
mtests_diss_tra_exp_II = mtests_diss_tra_exp_II[
  which(mtests_diss_tra_exp_II$`Pr(>|z|)`<0.05),]

# UNCERTAINTY * SIMILARITY * TRA
surv_df$res = interaction(surv_df$Uncertainty, surv_df$Diss, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_unc_diss_tra_res = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_unc_diss_tra = glht(model_unc_diss_tra_res, linfct=mcp(res="Tukey"))
summary(mcp_unc_diss_tra, test=adjusted(type="bonferroni"))
pq_unc_diss_tra = summary(mcp_unc_diss_tra, 
                           test=adjusted(type="bonferroni"))$test
mtests_unc_diss_tra_exp_II = cbind(
  pq_unc_diss_tra$coefficients, pq_unc_diss_tra$sigma, 
  pq_unc_diss_tra$tstat, pq_unc_diss_tra$pvalues)
error_unc_diss_tra = attr(pq_unc_diss_tra$pvalues, "error")
pname_unc_diss_tra = switch(
  mcp_unc_diss_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_unc_diss_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_unc_diss_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_unc_diss_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_unc_diss_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_unc_diss_tra$df==0,"z value","t value"), pname_unc_diss_tra)
mtests_unc_diss_tra_exp_II = as.data.frame(mtests_unc_diss_tra_exp_II)
mtests_unc_diss_tra_exp_II = mtests_unc_diss_tra_exp_II[
  which(mtests_unc_diss_tra_exp_II$`Pr(>|z|)`<0.05),]

# # PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS

file_name = "interaction_plot_sim_tra_unc_exp_II.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(model_unc_diss_tra, Diss ~ t_ra | Uncertainty, CIs = TRUE)
dev.off()

emm = emmeans(model_unc_diss_tra, specs=pairwise~t_ra|Uncertainty*Diss, 
              type="response", df=26)

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

# TABLES LATEX FOR SURVIE ~ UNCERTAINTY * SIMILARITY * TRA
xtable(anova(model_nul, model_unc_diss_tra))
xtable(anova(model_unc_diss_tra))
xtable(mtests_unc_exp_II)
xtable(mtests_diss_exp_II)
xtable(mtests_tra_exp_II)
xtable(mtests_unc_diss_exp_II)
xtable(mtests_unc_tra_exp_II)
xtable(mtests_diss_tra_exp_II)
xtable(mtests_unc_diss_tra_exp_II)
xtable(emm_df)

# model_unc_diss_tra = emmeans(model_unc_diss_tra, ~ Uncertainty * Diss * t_ra)
# contrast(model_unc_diss_tra, "eff", by = "Uncertainty")
# contrast(model_unc_diss_tra, "eff", by = "Diss")
# contrast(model_unc_diss_tra, "eff", by = "t_ra")
# contrast(model_unc_diss_tra, interaction = c("eff", "pairwise"),
#          adjust = multcomp_adjusted_type)
# xtable(contrast(model_unc_diss_tra, interaction = c("eff", "pairwise"),
#                 adjust = multcomp_adjusted_type))

################################################################################
# MODELE SURVIE ~ DENSITY * SIMILARITY * TRA
################################################################################

# SURVIVAL DATA ANALYZING
model_nul = coxph(survie ~ m_density * Diss * t_ra, data=surv_df)
model_den_diss_tra = coxph(survie ~ m_density * Diss * t_ra + 
                             frailty(Sujet, distribution="gamma"), data=surv_df)
anova(model_nul, model_den_diss_tra)
summary(model_den_diss_tra)
anova(model_den_diss_tra)

# MODELS DIAGNOSTICS
mres = resid(model_den_diss_tra, type = "martingale")
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

# # PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# # DENSITY
# mcp_den = glht(model_den_diss_tra, 
#                linfct=mcp(m_density="Tukey", 
#                           interaction_average = TRUE,
#                           covariate_average = TRUE))
# summary(mcp_den, test=adjusted(type=multcomp_adjusted_type))
# pq=summary(mcp_den)$test
# mtests_den_exp_II=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
# error=attr(pq$pvalues, "error")
# pname=switch(
#   mcp_den$alternativ,
#   less=paste("Pr(<", ifelse(mcp_den$df ==0, "z", "t"), ")", sep = ""),
#   greater=paste("Pr(>", ifelse(mcp_den$df==0, "z", "t"), ")", sep = ""),
#   two.sided=paste("Pr(>|", ifelse(mcp_den$df==0, "z", "t"), "|)",sep=""))
# colnames(mtests_den_exp_II) = c("Estimate", "Std. Error",
#                                 ifelse(mcp_den$df ==0, "z value", "t value"),
#                                 pname)

# # SIMILARITY
# mcp_diss = glht(model_den_diss_tra, 
#                 linfct=mcp(Diss="Tukey", 
#                            interaction_average = TRUE,
#                            covariate_average = TRUE))
# summary(mcp_diss, test=adjusted(type=multcomp_adjusted_type))
# pq = summary(mcp_diss)$test
# mtests_diss_exp_II = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
# error = attr(pq$pvalues, "error")
# pname = switch(
#   mcp_diss$alternativ,
#   less = paste("Pr(<", ifelse(mcp_diss$df ==0, "z", "t"), ")", sep = ""),
#   greater = paste("Pr(>", ifelse(mcp_diss$df == 0, "z", "t"), ")", sep=""),
#   two.sided=paste("Pr(>|", ifelse(mcp_diss$df==0, "z", "t"), "|)", sep=""))
# colnames(mtests_diss_exp_II) = c("Estimate", "Std. Error",
#                                  ifelse(mcp_diss$df==0, "z value", "t value"),
#                                  pname)

# # TRA
# mcp_tra = glht(model_den_diss_tra, 
#                linfct=mcp(t_ra="Tukey", 
#                           interaction_average = TRUE,
#                           covariate_average = TRUE))
# summary(mcp_tra,test=adjusted(type=multcomp_adjusted_type))
# pq = summary(mcp_tra)$test
# mtests_tra_exp_II = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
# error = attr(pq$pvalues, "error")
# pname = switch(
#   mcp_tra$alternativ,
#   less = paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""),
#   greater = paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep=""),
#   two.sided=paste("Pr(>|", ifelse(mcp_tra$df==0, "z", "t"), "|)", sep=""))
# colnames(mtests_tra_exp_II) = c("Estimate", "Std. Error",
#                                 ifelse(mcp_tra$df==0, "z value", "t value"),
#                                 pname)

# DENSITY * SIMILARITY
surv_df$res = interaction(surv_df$m_density, surv_df$Diss)
surv_df$res = factor(surv_df$res)
model_den_diss = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_den_diss = glht(model_den_diss, linfct=mcp(res="Tukey"))
summary(mcp_den_diss, test=adjusted(type="bonferroni"))
pq_den_diss = summary(mcp_den_diss, test=adjusted(type="bonferroni"))$test
mtests_den_diss_exp_II = cbind(pq_den_diss$coefficients, pq_den_diss$sigma, 
                               pq_den_diss$tstat, pq_den_diss$pvalues)
error_den_diss = attr(pq_den_diss$pvalues, "error")
pname_den_diss = switch(
  mcp_den_diss$alternativ,
  less = paste("Pr(<", ifelse(mcp_den_diss$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_den_diss$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_den_diss$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_den_diss_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_den_diss$df==0,"z value","t value"), pname_den_diss)
mtests_den_diss_exp_II = as.data.frame(mtests_den_diss_exp_II)
mtests_den_diss_exp_II = mtests_den_diss_exp_II[
  which(mtests_den_diss_exp_II$`Pr(>|z|)`<0.05),]

# DENSITY * TRA
surv_df$res = interaction(surv_df$m_density, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_den_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_den_tra = glht(model_den_tra, linfct=mcp(res="Tukey"))
summary(mcp_den_tra, test=adjusted(type="bonferroni"))
pq_den_tra = summary(mcp_den_tra, test=adjusted(type="bonferroni"))$test
mtests_den_tra_exp_II = cbind(pq_den_tra$coefficients, pq_den_tra$sigma, 
                              pq_den_tra$tstat, pq_den_tra$pvalues)
error_den_tra = attr(pq_den_tra$pvalues, "error")
pname_den_tra = switch(
  mcp_den_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_den_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_den_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_den_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_den_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_den_tra$df==0,"z value","t value"), pname_den_tra)
mtests_den_tra_exp_II = as.data.frame(mtests_den_tra_exp_II)
mtests_den_tra_exp_II = mtests_den_tra_exp_II[
  which(mtests_den_tra_exp_II$`Pr(>|z|)`<0.05),]

# SIMILARITY * TRA
surv_df$res = interaction(surv_df$Diss, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_diss_tra = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_diss_tra = glht(model_diss_tra, linfct=mcp(res="Tukey"))
summary(mcp_diss_tra, test=adjusted(type="bonferroni"))
pq_diss_tra = summary(mcp_diss_tra, test=adjusted(type="bonferroni"))$test
mtests_diss_tra_exp_II = cbind(pq_diss_tra$coefficients, pq_diss_tra$sigma, 
                               pq_diss_tra$tstat, pq_diss_tra$pvalues)
error_diss_tra = attr(pq_diss_tra$pvalues, "error")
pname_diss_tra = switch(
  mcp_diss_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_diss_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_diss_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_diss_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_diss_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_diss_tra$df==0,"z value","t value"), pname_diss_tra)
mtests_diss_tra_exp_II = as.data.frame(mtests_diss_tra_exp_II)
mtests_diss_tra_exp_II = mtests_diss_tra_exp_II[
  which(mtests_diss_tra_exp_II$`Pr(>|z|)`<0.05),]

# DENSITY * SIMILARITY * TRA
surv_df$res = interaction(surv_df$m_density, surv_df$Diss, surv_df$t_ra)
surv_df$res = factor(surv_df$res)
model_den_diss_tra_res = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_den_diss_tra = glht(model_den_diss_tra_res, linfct=mcp(res="Tukey"))
summary(mcp_den_diss_tra, test=adjusted(type="bonferroni"))
pq_den_diss_tra = summary(mcp_den_diss_tra, 
                          test=adjusted(type="bonferroni"))$test
mtests_den_diss_tra_exp_II = cbind(
  pq_den_diss_tra$coefficients, pq_den_diss_tra$sigma, 
  pq_den_diss_tra$tstat, pq_den_diss_tra$pvalues)
error_den_diss_tra = attr(pq_den_diss_tra$pvalues, "error")
pname_den_diss_tra = switch(
  mcp_den_diss_tra$alternativ,
  less = paste("Pr(<", ifelse(mcp_den_diss_tra$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(
    mcp_den_diss_tra$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(
    mcp_den_diss_tra$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_den_diss_tra_exp_II) = c(
  "Estimate", "Std. Error", ifelse(
    mcp_den_diss_tra$df==0,"z value","t value"), pname_den_diss_tra)
mtests_den_diss_tra_exp_II = as.data.frame(mtests_den_diss_tra_exp_II)
mtests_den_diss_tra_exp_II = mtests_den_diss_tra_exp_II[
  which(mtests_den_diss_tra_exp_II$`Pr(>|z|)`<0.05),]

# # PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS
# file_name = "interaction_plot_sim_tra_den_exp_II.pdf"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# pdf(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# emmip(model_den_diss_tra, Diss ~ t_ra | m_density, CIs = TRUE)
# dev.off()
# 
# emm = emmeans(model_den_diss_tra, specs=pairwise~t_ra|m_density*Diss, 
#               type="response", df=26)
# 
# emm_df = emm$contrasts %>%
#   summary(infer = TRUE) %>%
#   as.data.frame()

# TABLES LATEX FOR SURVIE ~ DENSITY * SIMILARITY * TRA
xtable(anova(model_nul, model_den_diss_tra))
xtable(anova(model_den_diss_tra))
# xtable(mtests_den_exp_II)
# xtable(mtests_diss_exp_II)
# xtable(mtests_tra_exp_II)
xtable(mtests_den_diss_exp_II)
xtable(mtests_den_tra_exp_II)
xtable(mtests_diss_tra_exp_II)
xtable(mtests_den_diss_tra_exp_II)
# xtable(emm_df)

# model_den_diss_tra = emmeans(model_den_diss_tra, ~ m_density * Diss * t_ra)
# contrast(model_den_diss_tra, "eff", by = "m_density")
# contrast(model_den_diss_tra, "eff", by = "Diss")
# contrast(model_den_diss_tra, "eff", by = "t_ra")
# contrast(model_den_diss_tra, interaction = c("eff", "pairwise"),
#          adjust = multcomp_adjusted_type)
# xtable(contrast(model_den_diss_tra, interaction = c("eff", "pairwise"),
#                 adjust = multcomp_adjusted_type))

################################################################################
# PLOTTING SURVIVAL FUNCTION, CUMULATIVE HAZARD FUNCTION AND CUMULATIVE 
# DISTRIBUTION FUNCTION
################################################################################

################################################################################
# CUMULATIVE DISTRIBUTION FUNCTION EN COULEUR ET TRAIT POINTILLE
################################################################################

# # survie ~ m_ppo
# c_mppo = survfit(survie ~ m_ppo, data = surv_df)
# file_name = "cdf_mppo_exp_II.eps"
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
# # survie ~ m_td
# c_mtd = survfit(survie ~ m_td, data = surv_df)
# file_name = "cdf_mtd_exp_II.eps"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# postscript(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(c_mtd, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
#      main = "", xlab = "Time (sec)", ylab = "F(t)", 
#      las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
# legend("bottomright", c("20ms", "60ms", "100ms"), lty=c(1,2,3), 
#        lwd=c(3,3,3), col = c("green", "red", "blue"), cex=1.6)
# dev.off()
# 
# # survie ~ t_ra
# c_tra = survfit(survie ~ t_ra, data = surv_df)
# file_name = "cdf_tra_exp_II.eps"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# postscript(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(c_tra, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
#      main = "", xlab = "Time (sec)", ylab = "F(t)", 
#      las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
# legend("bottomright", c("20ms", "60ms", "100ms"), lty=c(1,2,3), 
#        lwd=c(3,3,3), col = c("green", "red", "blue"), cex=1.6)
# dev.off()
# 
# # survie ~ Similarity
# c_diss = survfit(survie ~ Diss, data = surv_df)
# file_name = "cdf_similarity_exp_II.eps"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# postscript(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(c_diss, col = c(
#   "darkgreen","forestgreen","black", "green3", "green1"), 
#   lty=c(3,2,1,4,5), lwd=c(3,3,3),
#   main = "", xlab = "Time (sec)", ylab = "F(t)",
#   las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
# legend("bottomright", c("-80", "-40", "0", "40", "80"),
#        lty=c(3,2,1,4,5), lwd=c(3,3,3), col = c(
#          "darkgreen","forestgreen","black", "green3", "green1"), cex=1.6)
# dev.off()
# 
# # survie ~ Uncertainty
# c_uncertainty = survfit(survie ~ Uncertainty, data = surv_df)
# file_name = "cdf_uncertainty_exp_II.eps"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# postscript(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(c_uncertainty, col = c("green", "red", "blue"), 
#      lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", ylab = "F(t)",
#      las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
# legend("bottomright", c("29", "115", "463"), lty=c(1,2,3), lwd=c(3,3,3), 
#        col = c("green", "red", "blue"), cex=1.6)
# dev.off()
