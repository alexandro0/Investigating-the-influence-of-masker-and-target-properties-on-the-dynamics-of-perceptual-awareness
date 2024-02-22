################################################################################
# SCRIPT FOR ANALYZING SURVIVAL DATA FROM BEHAVIORAL EXPERIMENT DURING PHD 
# THESIS FOR EXPERIMENT I
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
article_target = "article_1_PA/PloSOne/maxi_dossier/Exp_I/dynamics_stats"
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
setwd("data/test/testA/sujet_marseille/")
surv_df=read.table("surv.csv",header=TRUE,sep=',',dec=".",fileEncoding="utf-8")

################################################################################
# DATA CLEANING/RESHAPING
################################################################################

surv_df$RT = surv_df$RT/1000
surv_df$m_ppo = factor(surv_df$m_ppo)
surv_df$m_td = factor(surv_df$m_td)
surv_df$t_td = factor(surv_df$t_td)
surv_df$Similarity = factor(surv_df$Similarity)
surv_df$m_density = factor(surv_df$m_density)
surv_df$Uncertainty = factor(surv_df$Uncertainty)
surv_df = surv_df[which(surv_df$Sujet != 2 & surv_df$Sujet != 6 &
                          surv_df$Sujet != 9),]
surv_df = surv_df[which(surv_df$Bloc != 1),] 
survie = Surv(surv_df$RT, surv_df$Surv)

# coln_original = c(
#   "X","Sujet","Bloc","Stim","Hits","FA","Miss","RC","RT","Similarity","m_ppo","m_td",
#   "m_iti","t_pi","t_ra","t_td","m_density","Uncertainty","Surv")
# coln_modifed = c(
#   "X","Sujet","Bloc","Stim","Hits","FA","Miss","RC","RT","Similarity","FPO",
#   "MTD","MITI","TPI","TRA","TTD","Density","Uncertainty","Surv")
# colnames(surv_df) = coln_modifed

################################################################################
# MODELE SURVIE ~ MPPO * MTD * TTD
################################################################################

# SURVIVAL DATA ANALYZING
model_nul = coxph(survie ~ m_ppo * m_td * t_td, data=surv_df)
model = coxph(survie ~ m_ppo * m_td * t_td + 
                frailty(Sujet, distribution="gamma"), data=surv_df)
anova(model_nul, model)
summary(model)
anova(model)

# MODELS DIAGNOSTICS
mres = resid(model, type = "martingale")
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
mcp_mppo = glht(model, linfct=mcp(m_ppo="Tukey", interaction_average = TRUE,
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
mcp_mtd = glht(model, linfct=mcp(m_td="Tukey", interaction_average = TRUE,
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
mcp_ttd = glht(model, linfct=mcp(t_td="Tukey", interaction_average = TRUE,
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

# MPPO * MTD
surv_df$res = interaction(surv_df$m_ppo, surv_df$m_td)
surv_df$res = factor(surv_df$res)
model_mppo_mtd = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
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
colnames(mtests_mppo_mtd_exp_I) = c("Estimate", "Std. Error", 
                               ifelse(mcp_mppo_mtd$df==0,"z value","t value"),
                               pname_mppo_mtd)
mtests_mppo_mtd_exp_I = as.data.frame(mtests_mppo_mtd_exp_I)
mtests_mppo_mtd_exp_I = mtests_mppo_mtd_exp_I[
  which(mtests_mppo_mtd_exp_I$`Pr(>|z|)`<0.05),]

# MPPO * TTD
surv_df$res = interaction(surv_df$m_ppo, surv_df$t_td)
surv_df$res = factor(surv_df$res)
model_mppo_ttd = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
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
mtests_mppo_ttd_exp_I = as.data.frame(mtests_mppo_ttd_exp_I)
mtests_mppo_ttd_exp_I = mtests_mppo_ttd_exp_I[
  which(mtests_mppo_ttd_exp_I$`Pr(>|z|)`<0.05),]

# MTD * TTD
surv_df$res = interaction(surv_df$m_td, surv_df$t_td)
surv_df$res = factor(surv_df$res)
model_mtd_ttd = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
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
mtests_mtd_ttd_exp_I = as.data.frame(mtests_mtd_ttd_exp_I)
mtests_mtd_ttd_exp_I = mtests_mtd_ttd_exp_I[
  which(mtests_mtd_ttd_exp_I$`Pr(>|z|)`<0.05),]

# MPPO * MTD * TTD
surv_df$res = interaction(surv_df$m_ppo, surv_df$m_td, surv_df$t_td)
surv_df$res = factor(surv_df$res)
model_mppo_mtd_ttd = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
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
mtests_mppo_mtd_ttd_exp_I = as.data.frame(mtests_mppo_mtd_ttd_exp_I)
mtests_mppo_mtd_ttd_exp_I = mtests_mppo_mtd_ttd_exp_I[
  which(mtests_mppo_mtd_ttd_exp_I$`Pr(>|z|)`<0.05),]

# TABLES LATEX FOR SURVIE ~ MPPO * MTD * TTD
xtable(anova(model_nul, model))
xtable(anova(model))
xtable(mtests_mppo_exp_I)
xtable(mtests_mtd_exp_I)
xtable(mtests_ttd_exp_I)
xtable(mtests_mppo_mtd_exp_I)
xtable(mtests_mppo_ttd_exp_I)
xtable(mtests_mtd_ttd_exp_I)
xtable(mtests_mppo_mtd_ttd_exp_I)

# model.emm = emmeans(model, ~ m_ppo * m_td * t_td)
# contrast(model.emm, "eff", by = "m_ppo")
# contrast(model.emm, "eff", by = "m_td")
# contrast(model.emm, "eff", by = "t_td")
# contrast(model.emm, interaction = c("eff", "pairwise"), 
#          adjust = multcomp_adjusted_type)
# xtable(contrast(model.emm, interaction = c("eff", "pairwise"), 
#                 adjust = multcomp_adjusted_type))

################################################################################
# MODELE SURVIE ~ MPPO * SIMILARITY
################################################################################

# SURVIVAL DATA ANALYZING
model_nul_mppo_sim = coxph(survie ~ m_ppo * Similarity, data=surv_df)
model_mppo_sim = coxph(survie ~ m_ppo * Similarity +
                     frailty(Sujet, distribution = "gamma"), data=surv_df)
anova(model_nul_mppo_sim, model_mppo_sim)
summary(model_mppo_sim)
anova(model_mppo_sim)

# MODELS DIAGNOSTICS
mres = resid(model_mppo_sim, type = "martingale")
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
plot(0,0,lty=1,type='n', xlim=c(0,3), ylim=c(0,3), main="",
     xlab="Cox-Snell residuals", ylab="Estimated Cumulative Hazard")
lines(rsurv$time, -log(rsurv$surv), type="s")
lines(c(0,3),c(0,3), col="red")

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# MPPO
mcp_mppo = glht(model_mppo_sim, linfct=mcp(
  m_ppo="Tukey", interaction_average = TRUE, covariate_average = TRUE))
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

# SIMILARITY
mcp_sim = glht(model_mppo_sim, linfct=mcp(
  Similarity="Tukey", interaction_average = TRUE, covariate_average = TRUE))
summary(mcp_sim, test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_sim)$test
mtests_sim_exp_I = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_sim$alternativ,
  less = paste("Pr(<", ifelse(mcp_sim$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_sim$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_sim$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_sim_exp_I) = c("Estimate", "Std. Error",
                               ifelse(mcp_sim$df==0, "z value", "t value"),
                               pname)

# MPPO * SIMILARITY
surv_df$res = interaction(surv_df$m_ppo, surv_df$Similarity)
surv_df$res = factor(surv_df$res)
model_mppo_sim = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_mppo_sim = glht(model_mppo_sim, linfct=mcp(res="Tukey"))
summary(mcp_mppo_sim, test=adjusted(type="bonferroni"))
pq_mppo_sim = summary(mcp_mppo_sim, test=adjusted(type="bonferroni"))$test
mtests_mppo_sim_exp_I = cbind(pq_mppo_sim$coefficients, pq_mppo_sim$sigma, 
                              pq_mppo_sim$tstat, pq_mppo_sim$pvalues)
error_mppo_sim = attr(pq_mppo_sim$pvalues, "error")
pname_mppo_sim = switch(
  mcp_mppo_sim$alternativ,
  less = paste("Pr(<", ifelse(mcp_mppo_sim$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_mppo_sim$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_mppo_sim$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_mppo_sim_exp_I) = c("Estimate", "Std. Error", 
                                    ifelse(mcp_mppo_sim$df==0,"z value","t value"),
                                    pname_mppo_sim)
mtests_mppo_sim_exp_I = as.data.frame(mtests_mppo_sim_exp_I)
mtests_mppo_sim_exp_I = mtests_mppo_sim_exp_I[
  which(mtests_mppo_sim_exp_I$`Pr(>|z|)`<0.05),]

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS
model_mppo_sim = coxph(survie ~ m_ppo * Similarity +
                         frailty(Sujet, distribution = "gamma"), data=surv_df)

file_name = "interaction_plot_mfpo_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(model_mppo_sim, m_ppo ~ Similarity)
dev.off()

emm = emmeans(model_mppo_sim, specs=pairwise~m_ppo|Similarity, type="response", df=14)

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_mfpo_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(model_mppo_sim, ~m_ppo|Similarity, type="response")

file_name = "pwpp_mfpo_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

# SURVIE ~ MPPO * SIMILARITY
xtable(anova(model_nul_mppo_sim, model_mppo_sim))
xtable(anova(model_mppo_sim))
xtable(mtests_mppo_exp_I)
xtable(mtests_sim_exp_I)
xtable(mtests_mppo_sim_exp_I)
xtable(emm_df)

# model_Similarity.emm = emmeans(model_Similarity, ~ m_ppo * Similarity)
# contrast(model_Similarity.emm, "eff", by = "m_ppo")
# contrast(model_Similarity.emm, "eff", by = "Similarity")
# contrast(model_Similarity.emm, interaction = c("eff", "pairwise"), 
#          adjust = multcomp_adjusted_type, df=14)
# xtable(contrast(model_Similarity.emm, interaction = c("eff", "pairwise"), 
#                 adjust = multcomp_adjusted_type))

################################################################################
# MODELE SURVIE ~ UNCERTAINTY * SIMILARITY
################################################################################

# SURVIVAL DATA ANALYZING
model_nul_unc_sim = coxph(survie ~ Uncertainty * Similarity, data=surv_df)
model_unc_sim = coxph(survie ~ Uncertainty * Similarity + 
                        frailty(Sujet, distribution = "gamma"), data=surv_df)
anova(model_nul_unc_sim, model_unc_sim)
summary(model_unc_sim)
anova(model_unc_sim)

# MODELS DIAGNOSTICS
mres = resid(model_unc_sim, type = "martingale")
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
plot(0,0,lty=1,type='n', xlim=c(0,12), ylim=c(0,12), main="",
     xlab="Cox-Snell residuals", ylab="Estimated Cumulative Hazard")
lines(rsurv$time, -log(rsurv$surv), type="s")
lines(c(0,9),c(0,9), col="red")

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# UNCERTAINTY
mcp_unc = glht(model_unc_sim, linfct=mcp(
  Uncertainty="Tukey", interaction_average = TRUE, covariate_average = TRUE))
summary(mcp_unc, test=adjusted(type=multcomp_adjusted_type))
pq=summary(mcp_unc)$test
mtests_unc_exp_I=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
  mcp_unc$alternativ,
  less=paste("Pr(<", ifelse(mcp_unc$df ==0, "z", "t"), ")", sep = ""),
  greater=paste("Pr(>", ifelse(mcp_unc$df==0, "z", "t"), ")", sep = ""),
  two.sided=paste("Pr(>|", ifelse(mcp_unc$df==0, "z", "t"), "|)",sep=""))
colnames(mtests_unc_exp_I) = c("Estimate", "Std. Error",
                                ifelse(mcp_unc$df ==0, "z value", "t value"),
                                pname)

# SIMILARITY
mcp_sim = glht(model_unc_sim, linfct=mcp(
  Similarity="Tukey", interaction_average = TRUE, covariate_average = TRUE))
summary(mcp_sim, test=adjusted(type=multcomp_adjusted_type))
pq = summary(mcp_sim)$test
mtests_sim_exp_I = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error = attr(pq$pvalues, "error")
pname = switch(
  mcp_sim$alternativ,
  less = paste("Pr(<", ifelse(mcp_sim$df ==0, "z", "t"), ")", sep = ""),
  greater = paste("Pr(>", ifelse(mcp_sim$df == 0, "z", "t"), ")", sep=""),
  two.sided=paste("Pr(>|", ifelse(mcp_sim$df==0, "z", "t"), "|)", sep=""))
colnames(mtests_sim_exp_I) = c("Estimate", "Std. Error",
                               ifelse(mcp_sim$df==0, "z value", "t value"),
                               pname)

# UNCERTAINTY * SIMILARITY
surv_df$res = interaction(surv_df$Uncertainty, surv_df$Similarity)
surv_df$res = factor(surv_df$res)
model_unc_sim = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_unc_sim = glht(model_unc_sim, linfct=mcp(res="Tukey"))
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
colnames(mtests_unc_sim_exp_I) = c("Estimate", "Std. Error", 
                                    ifelse(mcp_unc_sim$df==0,"z value","t value"),
                                    pname_unc_sim)
mtests_unc_sim_exp_I = as.data.frame(mtests_unc_sim_exp_I)
mtests_unc_sim_exp_I = mtests_unc_sim_exp_I[
  which(mtests_unc_sim_exp_I$`Pr(>|z|)`<0.05),]

# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS
model_unc_sim = coxph(survie ~ Uncertainty * Similarity + 
                        frailty(Sujet, distribution = "gamma"), data=surv_df)

file_name = "interaction_plot_uncertainty_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
emmip(model_unc_sim, Uncertainty ~ Similarity)
dev.off()

emm = emmeans(model_unc_sim, specs=pairwise~Uncertainty|Similarity, type="response", df=14)

emm_df = emm$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

file_name = "emmeans_uncertainty_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(emm, comparisons = TRUE)
dev.off()

emmn = emmeans(model_unc_sim, ~Uncertainty|Similarity, type="response")

file_name = "pwpp_uncertainty_similarity_exp_I.pdf"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
pdf(file = figure_to_save, width=14, height=6)
par(mar=c(4,5,0,0)+0.2)
pwpp(emmn)
dev.off()

# TABLES LATEX SURVIE ~ UNCERTAINTY * SIMILARITY
xtable(anova(model_nul_unc_sim, model_unc_sim))
xtable(anova(model_unc_sim))
xtable(mtests_unc_exp_I)
xtable(mtests_sim_exp_I)
xtable(mtests_unc_sim_exp_I)
xtable(emm_df)

# model_unc_sim = emmeans(model_unc_sim, ~ Uncertainty * Similarity, df=14)
# contrast(model_unc_sim, "eff", by = "Uncertainty")
# contrast(model_unc_sim, "eff", by = "Similarity")
# contrast(model_unc_sim, interaction = c("eff", "pairwise"), 
#          adjust = multcomp_adjusted_type, df=14)
# xtable(contrast(model_unc_sim, interaction = c("eff", "pairwise"),
#                 adjust = multcomp_adjusted_type))

################################################################################
# MODELE SURVIE ~ DENSITY * SIMILARITY
################################################################################

# SURVIVAL DATA ANALYZING
model_nul_den_sim = coxph(survie ~ m_density * Similarity, data=surv_df)
model_den_sim = coxph(survie ~ m_density * Similarity + 
                        frailty(Sujet, distribution = "gamma"), data=surv_df)
anova(model_nul_den_sim, model_den_sim)
summary(model_den_sim)
anova(model_den_sim)

# MODELS DIAGNOSTICS
mres = resid(model_den_sim, type = "martingale")
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
plot(0,0,lty=1,type='n', xlim=c(0,3), ylim=c(0,3), main="",
     xlab="Cox-Snell residuals", ylab="Estimated Cumulative Hazard")
lines(rsurv$time, -log(rsurv$surv), type="s")
lines(c(0,3),c(0,3), col="red")

# # PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH MULTCOMP
# # DENSITY
# mcp_den = glht(model_den_sim, linfct=mcp(
#   m_density="Tukey", interaction_average = TRUE, covariate_average = TRUE))
# summary(mcp_den, test=adjusted(type=multcomp_adjusted_type))
# pq=summary(mcp_den)$test
# mtests_den_exp_I=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
# error=attr(pq$pvalues, "error")
# pname=switch(
#   mcp_den$alternativ,
#   less=paste("Pr(<", ifelse(mcp_den$df ==0, "z", "t"), ")", sep = ""),
#   greater=paste("Pr(>", ifelse(mcp_den$df==0, "z", "t"), ")", sep = ""),
#   two.sided=paste("Pr(>|", ifelse(mcp_den$df==0, "z", "t"), "|)",sep=""))
# colnames(mtests_den_exp_I) = c("Estimate", "Std. Error",
#                                ifelse(mcp_den$df ==0, "z value", "t value"),
#                                pname)
# 
# # SIMILARITY
# mcp_sim = glht(model_den_sim, linfct=mcp(
#   Similarity="Tukey", interaction_average = TRUE, covariate_average = TRUE))
# summary(mcp_sim, test=adjusted(type=multcomp_adjusted_type))
# pq = summary(mcp_sim)$test
# mtests_sim_exp_I = cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
# error = attr(pq$pvalues, "error")
# pname = switch(
#   mcp_sim$alternativ,
#   less = paste("Pr(<", ifelse(mcp_sim$df ==0, "z", "t"), ")", sep = ""),
#   greater = paste("Pr(>", ifelse(mcp_sim$df == 0, "z", "t"), ")", sep=""),
#   two.sided=paste("Pr(>|", ifelse(mcp_sim$df==0, "z", "t"), "|)", sep=""))
# colnames(mtests_sim_exp_I) = c("Estimate", "Std. Error",
#                                ifelse(mcp_sim$df==0, "z value", "t value"),
#                                pname)

# DENSITY * SIMILARITY
surv_df$res = interaction(surv_df$m_density, surv_df$Similarity)
surv_df$res = factor(surv_df$res)
model_den_sim = coxph(
  survie ~ res + frailty(Sujet, distribution="gamma"), data=surv_df)
mcp_den_sim = glht(model_den_sim, linfct=mcp(res="Tukey"))
summary(mcp_den_sim, test=adjusted(type="bonferroni"))
pq_den_sim = summary(mcp_den_sim, test=adjusted(type="bonferroni"))$test
mtests_den_sim_exp_I = cbind(pq_den_sim$coefficients, pq_den_sim$sigma, 
                             pq_den_sim$tstat, pq_den_sim$pvalues)
error_den_sim = attr(pq_den_sim$pvalues, "error")
pname_den_sim = switch(
  mcp_den_sim$alternativ,
  less = paste("Pr(<", ifelse(mcp_den_sim$df ==0, "z", "t"), ")",sep=""),
  greater = paste("Pr(>", ifelse(mcp_den_sim$df == 0, "z", "t"), ")",sep=""),
  two.sided = paste("Pr(>|", ifelse(mcp_den_sim$df == 0,"z","t"),"|)",sep=""))
colnames(mtests_den_sim_exp_I) = c("Estimate", "Std. Error", 
                                   ifelse(mcp_den_sim$df==0,"z value","t value"),
                                   pname_den_sim)
mtests_den_sim_exp_I = as.data.frame(mtests_den_sim_exp_I)
mtests_den_sim_exp_I = mtests_den_sim_exp_I[
  which(mtests_den_sim_exp_I$`Pr(>|z|)`<0.05),]

# # PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH EMMEANS
# model_den_sim = coxph(survie ~ m_density * Similarity +
#                         frailty(Sujet, distribution = "gamma"), data=surv_df)
# 
# file_name = "interaction_plot_density_similarity_exp_I.pdf"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# pdf(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# emmip(model_den_sim, m_density ~ Similarity)
# dev.off()
# 
# emm = emmeans(model_den_sim, specs=pairwise~m_density|Similarity, type="response", df=14)
# 
# emm_df = emm$contrasts %>%
#   summary(infer = TRUE) %>%
#   as.data.frame()
# 
# file_name = "emmeans_density_similarity_exp_I.pdf"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# pdf(file = figure_to_save, width=7, height=5)
# par(mar=c(4,5,0,0)+0.2)
# plot(emm, comparisons = TRUE)
# dev.off()
# 
# emmn = emmeans(model_den_sim, ~m_density|Similarity, type="response")
# 
# file_name = "pwpp_density_similarity_exp_I.pdf"
# figure_to_save = str_c(article_path, article_target, file_name, sep="/")
# pdf(file = figure_to_save, width=14, height=6)
# par(mar=c(4,5,0,0)+0.2)
# pwpp(emmn)
# dev.off()

# TABLES LATEX SURVIE ~ UNCERTAINTY * SIMILARITY
model_den_sim = coxph(survie ~ m_density * Similarity + 
                        frailty(Sujet, distribution = "gamma"), data=surv_df)

xtable(anova(model_nul_den_sim, model_den_sim))
xtable(anova(model_den_sim))
# xtable(mtests_den_exp_I)
# xtable(mtests_sim_exp_I)
xtable(mtests_den_sim_exp_I)
# xtable(emm_df)

################################################################################
# PLOTTING SURVIVAL FUNCTION, CUMULATIVE HAZARD FUNCTION AND CUMULATIVE 
# DISTRIBUTION FUNCTION
################################################################################

################################################################################
# CUMULATIVE DISTRIBUTION FUNCTION EN COULEUR ET TRAIT POINTILLE
################################################################################

# survie ~ m_ppo
c_mppo = survfit(survie ~ m_ppo, data = surv_df)
file_name = "cdf_mppo_exp_I.eps"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
postscript(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(c_mppo, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),
     main = "", xlab = "Time (sec)", ylab = "F(t)",
     las=1, cex.axis=1.6, cex.lab=1.6, fun="hazard")
legend("bottomright", c("4fpo", "32fpo", "64fpo"), lty=c(1,2,3), lwd=c(3,3,3),
       col = c("green", "red", "blue"), cex=1.6)
dev.off()

# survie ~ m_td
c_mtd = survfit(survie ~ m_td, data = surv_df)
file_name = "cdf_mtd_exp_I.eps"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
postscript(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(c_mtd, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),
     main = "", xlab = "Time (sec)", ylab = "F(t)",
     las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
legend("bottomright", c("20ms", "60ms", "100ms"), lty=c(1,2,3),
       lwd=c(3,3,3), col = c("green", "red", "blue"), cex=1.6)
dev.off()

# survie ~ t_td
c_ttd = survfit(survie ~ t_td, data = surv_df)
file_name = "cdf_ttd_exp_I.eps"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
postscript(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(c_ttd, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),
     main = "", xlab = "Time (sec)", ylab = "F(t)",
     las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
legend("bottomright", c("20ms", "60ms", "100ms"), lty=c(1,2,3),
       lwd=c(3,3,3), col = c("green", "red", "blue"), cex=1.6)
dev.off()

# survie ~ Similarity
c_Similarity = survfit(survie ~ Similarity, data = surv_df)
file_name = "cdf_similarity_exp_I.eps"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
postscript(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(c_Similarity, col = c(
  "darkgreen","forestgreen","black", "green3", "green1"),
  lty=c(3,2,1,4,5), lwd=c(3,3,3),
  main = "", xlab = "Time (sec)", ylab = "F(t)",
  las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
legend("bottomright", c("-80", "-40", "0", "40", "80"),
       lty=c(3,2,1,4,5), lwd=c(3,3,3), col = c(
         "darkgreen","forestgreen","black", "green3", "green1"), cex=1.6)
dev.off()

# survie ~ Uncertainty
c_uncertainty = survfit(survie ~ Uncertainty, data = surv_df)
file_name = "cdf_uncertainty_exp_I.eps"
figure_to_save = str_c(article_path, article_target, file_name, sep="/")
postscript(file = figure_to_save, width=7, height=5)
par(mar=c(4,5,0,0)+0.2)
plot(c_uncertainty, col = c("green", "red", "blue"),
     lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", ylab = "F(t)",
     las=1, cex.axis=1.6, cex.lab=1.6, fun="F")
legend("bottomright", c("29", "115", "463"), lty=c(1,2,3), lwd=c(3,3,3),
       col = c("green", "red", "blue"), cex=1.6)
dev.off()
