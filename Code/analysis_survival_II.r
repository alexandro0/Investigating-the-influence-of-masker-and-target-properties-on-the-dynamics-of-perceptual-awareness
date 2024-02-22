################################################################################
# SCRIPT FOR ANALYZING SURVIVAL DATA FROM BEHAVIORAL EXPERIMENT DURING PHD 
# THESIS FOR THE SECOND EXPERIMENT (EXP II)
################################################################################

# CTRL + ALT + R

par(mar=c(4,4,0,0)+0.2)

################################################################################
# LIBRARY LOAD
################################################################################

library(Matrix)
library(MASS)
library(bdsmatrix)
library(lattice)
library(mvtnorm)
library(TH.data)
library(lme4)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(multcomp)
library(survival)
library(coxme)
library(frailtyEM)
library(xtable)

################################################################################
# DATA LOAD
################################################################################

setwd("/home/link/Documents/thèse_onera/experimentation/psychophysique/")
setwd("data/test/testB/sujet_marseille/")
data=read.table("data.csv",header=TRUE,sep=',',dec=".",fileEncoding="utf-16")
surv_df=read.table("surv.csv",header=TRUE,sep=',',dec=".",fileEncoding="utf-16")

################################################################################
# DATA CLEANING/RESHAPING
################################################################################

data$Stim = NULL 
data$X = NULL
# data$m_density = round(data$m_density)
tmp = data[which(data$RT != 0 & data$RT > 700 & data$Hits == 1),]
tpm = data[which(data$RT != 0 & data$RT > 700 & data$Hits == 1),]
tmp = tmp[which(tmp$Sujet!=1 & tmp$Sujet!=2 & tmp$Sujet!=6 & tmp$Sujet!=8),]
tmp = tmp[which(tmp$Bloc != 1),] 

surv_df$RT = surv_df$RT/1000
surv_df$m_ppo = factor(surv_df$m_ppo)
surv_df$m_td = factor(surv_df$m_td)
surv_df$t_ra = factor(surv_df$t_ra)
surv_df$Diss = factor(surv_df$Diss)
surv_df$m_density = factor(surv_df$m_density)
surv_df = surv_df[which(surv_df$Sujet != 1 & surv_df$Sujet != 2 
                        & surv_df$Sujet != 6 & surv_df$Sujet != 8),]
surv_df = surv_df[which(surv_df$Bloc != 1),]
survie = Surv(surv_df$RT, surv_df$Surv)

################################################################################
# CLASSICAL DATA ANALYZING
################################################################################

par(mfrow=c(2,2))
hist(tmp[,"RT"], main="", breaks=30, col=c("skyblue"), prob=TRUE, xlab="DT(ms)")
lines(density(tmp[,"RT"]))
abline(v=mean(tmp[,"RT"]), col="green")
gamm = rgamma(900, 1, 0.6)
hist(gamm, breaks=30, col=c("skyblue"), prob=TRUE, xlab="Gamma Distr.", main="")
lines(density(gamm))

boxplot(RT ~ Sujet, col=c("skyblue"), data = tpm, main="", xlab="Sujet", 
        ylab="Detection Time (ms)")
abline(h=mean(tpm[,"RT"]), col="green")

boxplot(RT~Bloc, col=c("skyblue"), data = tmp, xlab="Bloc", ylab="DT", main="")
tmp$Bloc = factor(tmp$Bloc)
bloc_reg = lm(RT~Bloc, data=tmp)
summary(bloc_reg)
anova(bloc_reg)
summary(glht(bloc_reg, linfct=mcp(Bloc="Tukey", interaction_average=TRUE)))

boxplot(RT ~ m_density, col=c("skyblue"), data = tmp, xlab="Mask Density", 
        ylab="DT", main="")
tmp$m_density = factor(tmp$m_density)
mden_reg = lm(RT ~ m_density, data=tmp)
summary(mden_reg)
anova(mden_reg)
summary(glht(mden_reg, linfct=mcp(m_density="Tukey", interaction_average=TRUE)))

boxplot(RT~m_ppo, col=c("skyblue"), data = tmp, xlab="MPPO", ylab="DT", main="")
tmp$m_ppo = factor(tmp$m_ppo)
mppo_reg = lm(RT ~ m_ppo, data=tmp)
summary(mppo_reg)
anova(mppo_reg)
summary(glht(mppo_reg, linfct=mcp(m_ppo="Tukey", interaction_average=TRUE)))

boxplot(RT ~ m_td, col=c("skyblue"), data = tmp, xlab="MTD", ylab="DT", main="")
tmp$m_td = factor(tmp$m_td)
mtd_reg = lm(RT ~ m_td, data=tmp)
summary(mtd_reg)
anova(mtd_reg)
summary(glht(mtd_reg, linfct=mcp(m_td="Tukey", interaction_average=TRUE)))

boxplot(RT ~ t_ra, col=c("skyblue"), data = tmp, xlab="TRA", ylab="DT", main="")
tmp$t_ra = factor(tmp$t_ra)
tra_reg = lm(RT ~ t_ra, data=tmp)
summary(tra_reg)
anova(tra_reg)
summary(glht(tra_reg, linfct=mcp(t_ra="Tukey", interaction_average=TRUE)))

boxplot(RT ~ Diss, names=c("S", "LD", "HD"), col=c("skyblue"), 
        xlab="Similarity", ylab="DT", data = tmp, main="")
tmp$Diss = factor(tmp$Diss)
diss_reg = lm(RT ~ Diss, data=tmp)
summary(diss_reg)
anova(diss_reg)
summary(glht(diss_reg, linfct=mcp(Diss="Tukey", interaction_average=TRUE)))

################################################################################
# SURVIVAL DATA ANALYZING
################################################################################

model_nul = coxph(survie ~ m_ppo * m_td * t_ra, data = surv_df)
model = coxph(survie ~ m_ppo * m_td * t_ra + 
                      frailty(Sujet, distribution = "gamma"), data=surv_df)
model_diss = coxph(survie ~ m_ppo * Diss * t_ra + 
                           frailty(Sujet, distribution = "gamma"), data=surv_df)

anova(model_nul, model)
summary(model)
anova(model)

anova(model_nul, model_diss)
summary(model_diss)
anova(model_diss)

################################################################################
# MODELS DIAGNOSTICS
################################################################################

mres = resid(model, type = "martingale")
csres = surv_df$Surv - mres
exp_dist = rexp(1000, rate=1)
par(mfrow=c(2,2))
hist(csres, breaks = 30, col=c("lightblue"), prob=TRUE, main="", 
     xlab="Cox-Snell Residuals", ylab="Density", ylim=c(0,1), xlim=c(0,6))
lines(density(csres), col="red")
hist(exp_dist, breaks = 30, col=c("lightblue"), prob=TRUE, main="", 
     xlab="Exponential Distribution", ylab="Density", ylim=c(0,1), xlim=c(0,6))
lines(density(exp_dist), col="red")
rsurv = survfit(Surv(csres, surv_df$Surv) ~ 1, type="fleming-harrington")
plot(0, 0, lty=1, type='n', xlim=c(0,3), ylim=c(0,3), 
     main="Diagnostic plot of the frailty model in Experiment II", 
     xlab="Cox-Snell residuals", 
     ylab="Estimated Cumulative Hazard")
lines(rsurv$time, -log(rsurv$surv), type="s")
lines(c(0,3),c(0,3),col="red")

mres = resid(model_diss, type = "martingale")
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

################################################################################
# GENERATION OF LATEX DATA TABLE
################################################################################

xtable(anova(model_nul, model))
xtable(anova(model))

xtable(anova(model_nul, model_diss))
xtable(anova(model_diss))

################################################################################
# PERFORMING MULTIPLE POST-HOCS COMPARISONS WITH BONFERO=RONI CORRECTION
################################################################################

# summary(glht(model, linfct=mcp(m_ppo="Tukey", interaction_average=TRUE)))
# summary(glht(model, linfct=mcp(m_td="Tukey", interaction_average=TRUE)))
# summary(glht(model, linfct=mcp(t_ra="Tukey", interaction_average=TRUE)))
# summary(glht(model_diss, linfct=mcp(Diss="Tukey", interaction_average=TRUE)))

mcp_mppo = glht(model, linfct=mcp(m_ppo="Tukey", interaction_average = TRUE, 
                                  covariate_average = TRUE))
summary(mcp_mppo,test=adjusted(type="bonferroni"))
pq=summary(mcp_mppo)$test
mtests=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
        mcp_mppo$alternativ, 
        less=paste("Pr(<", ifelse(mcp_mppo$df ==0, "z", "t"), ")", sep = ""),
        greater=paste("Pr(>", ifelse(mcp_mppo$df==0, "z", "t"), ")", sep = ""),
        two.sided=paste("Pr(>|", ifelse(mcp_mppo$df==0, "z", "t"), "|)",sep=""))
colnames(mtests) = c("Estimate", "Std. Error", 
                     ifelse(mcp_mppo$df ==0, "z value", "t value"), pname)
xtable(mtests)

mcp_mtd = glht(model, linfct=mcp(m_td="Tukey", interaction_average = TRUE, 
                                 covariate_average = TRUE))
summary(mcp_mtd,test=adjusted(type="bonferroni"))
pq=summary(mcp_mtd)$test
mtests=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
        mcp_mtd$alternativ, 
        less=paste("Pr(<", ifelse(mcp_mtd$df ==0, "z", "t"), ")", sep = ""),
        greater=paste("Pr(>", ifelse(mcp_mtd$df == 0, "z", "t"), ")", sep=""),
        two.sided=paste("Pr(>|",ifelse(mcp_mtd$df==0, "z", "t"), "|)", sep=""))
colnames(mtests) = c("Estimate", "Std. Error", 
                     ifelse(mcp_mtd$df ==0, "z value", "t value"), pname)
xtable(mtests)

mcp_tra = glht(model, linfct=mcp(t_ra="Tukey", interaction_average = TRUE, 
                                 covariate_average = TRUE))
summary(mcp_tra,test=adjusted(type="bonferroni"))
pq=summary(mcp_tra)$test
mtests=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
        mcp_tra$alternativ, 
        less=paste("Pr(<", ifelse(mcp_tra$df ==0, "z", "t"), ")", sep = ""),
        greater=paste("Pr(>", ifelse(mcp_tra$df == 0, "z", "t"), ")", sep = ""),
        two.sided=paste("Pr(>|", ifelse(mcp_tra$df==0, "z", "t"), "|)", sep=""))
colnames(mtests) = c("Estimate", "Std. Error", 
                     ifelse(mcp_tra$df ==0, "z value", "t value"), pname)
xtable(mtests)

mcp_diss = glht(model_diss, linfct=mcp(Diss="Tukey", interaction_average = TRUE, 
                                       covariate_average = TRUE))
summary(mcp_diss,test=adjusted(type="bonferroni"))
pq=summary(mcp_diss)$test
mtests=cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error=attr(pq$pvalues, "error")
pname=switch(
        mcp_diss$alternativ, 
        less=paste("Pr(<", ifelse(mcp_diss$df ==0, "z", "t"), ")", sep = ""),
        greater=paste("Pr(>", ifelse(mcp_diss$df==0, "z", "t"), ")", sep = ""),
        two.sided=paste("Pr(>|", ifelse(mcp_diss$df==0, "z", "t"), "|)",sep=""))
colnames(mtests) = c("Estimate", "Std. Error", 
                     ifelse(mcp_diss$df ==0, "z value", "t value"), pname)
xtable(mtests)

################################################################################
# PLOTTING SURVIVAL FUNCTION, CUMULATIVE HAZARD FUNCTION AND CUMULATIVE 
# DISTRIBUTION FUNCTION
################################################################################

# par(mar=c(0,0,0,0)+0.2)
# par(mar=c(4,4,0,0)+0.2)
# par(mfrow=c(2,2))
# pch=c(16,15,14)

# library(stringr)
# paper_path_1 = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/"
# paper_path_2 = "v21/figures/figure_courbes_frailtymodel_all/"
# figure_to_save = str_c(paper_path_1, paper_path_2)
# name_fig = "SF_EXP_II_MASKER_SPECTRAL_DENSITY_NB.jpeg"
# to_save = str_c(figure_to_save, name_fig)

################################################################################
# SURVIVAL FUNCTION EN NOIR ET BLANC ET TRAIT POINTILLE
################################################################################
c_mppo = survfit(survie ~ m_ppo, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/SF_EXP_II_MASKER_SPECTRAL_DENSITY_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mppo, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Survival Function (%)", las=1, cex.axis=1.3, cex.lab=1.3)
legend("topright", c("Low", "Medium", "High"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_mtd = survfit(survie ~ m_td, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/SF_EXP_II_MASKER_TONE_DURATION_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mtd, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Survival Function (%)", las=1, cex.axis=1.3, cex.lab=1.3)
legend("topright", c("Short","Intermediate","Long"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_tra = survfit(survie ~ t_ra, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/SF_EXP_II_TARGET_RATE_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_tra, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Survival Function (%)", las=1, cex.axis=1.3, cex.lab=1.3)
legend("topright", c("5Hz","10Hz","20Hz"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_diss = survfit(survie ~ Diss, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/SF_EXP_II_MASKER_TARGET_SIMILARITY_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_diss, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Survival Function (%)", las=1, cex.axis=1.3, cex.lab=1.3)
legend("topright", c("High Similarity", "Medium Similarity", "Low Similarity"), 
       lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

################################################################################
# SURVIVAL FUNCTION EN COULEUR ET TRAIT POINTILLE
################################################################################
c_mppo = survfit(survie ~ m_ppo, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/SF_EXP_II_MASKER_SPECTRAL_DENSITY_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mppo, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3), 
     main = "", xlab = "Time (sec)", ylab = "Survival Function (%)", las=1, 
     cex.axis=1.3, cex.lab=1.3)
legend("topright", c("Low", "Medium", "High"), lty=c(1,2,3), lwd=c(3,3,3), 
       col = c("green", "red", "blue"))
dev.off()

c_mtd = survfit(survie ~ m_td, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/SF_EXP_II_MASKER_TONE_DURATION_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mtd, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3), 
     main = "", xlab = "Time (sec)", ylab = "Survival Function (%)", las=1, 
     cex.axis=1.3, cex.lab=1.3)
legend("topright", c("Short", "Intermediate", "Long"), lty=c(1,2,3), 
       lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

c_tra = survfit(survie ~ t_ra, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/SF_EXP_II_TARGET_RATE_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_tra, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
     main = "", xlab = "Time (sec)", ylab = "Survival Function (%)", las=1, 
     cex.axis=1.3, cex.lab=1.3)
legend("topright", c("5Hz","10Hz","20Hz"), lty=c(1,2,3), 
       lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

c_diss = survfit(survie ~ Diss, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/SF_EXP_II_MASKER_TARGET_SIMILARITY_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_diss, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3), 
     main = "", xlab = "Time (sec)", ylab = "Survival Function (%)", las=1, 
     cex.axis=1.3, cex.lab=1.3)
legend("topright", c("High Similarity", "Medium Similarity", "Low Similarity"), 
       lty=c(1,2,3), lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

################################################################################
# CUMULATIVE HAZARD FUNCTION EN NOIR ET BLANC ET TRAIT POINTILLE
################################################################################
c_mppo = survfit(survie ~ m_ppo, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CHF_EXP_II_MASKER_SPECTRAL_DENSITY_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mppo, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Cumulative Hazard Function", las=1, cex.axis=1.3, cex.lab=1.3, 
     fun="cumhaz")
legend("topleft", c("Low", "Medium", "High"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_mtd = survfit(survie ~ m_td, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CHF_EXP_II_MASKER_TONE_DURATION_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mtd, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Cumulative Hazard Function", las=1, cex.axis=1.3, cex.lab=1.3, 
     fun="cumhaz")
legend("topleft", c("Short","Intermediate","Long"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_tra = survfit(survie ~ t_ra, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CHF_EXP_II_TARGET_RATE_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_tra, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Cumulative Hazard Function", las=1, cex.axis=1.3, cex.lab=1.3, 
     fun="cumhaz")
legend("topleft", c("5Hz","10Hz","20Hz"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_diss = survfit(survie ~ Diss, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CHF_EXP_II_MASKER_TARGET_SIMILARITY_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_diss, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Cumulative Hazard Function", las=1, cex.axis=1.3, cex.lab=1.3, 
     fun="cumhaz")
legend("topleft", c("High Similarity", "Medium Similarity", "Low Similarity"), 
       lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

################################################################################
# CUMULATIVE HAZARD FUNCTION EN COULEUR ET TRAIT POINTILLE
################################################################################
c_mppo = survfit(survie ~ m_ppo, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CHF_EXP_II_MASKER_SPECTRAL_DENSITY_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mppo, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3), 
     main = "", xlab = "Time (sec)", ylab = "Cumulative Hazard Function", las=1, 
     cex.axis=1.3, cex.lab=1.3, fun="cumhaz")
legend("topleft", c("Low", "Medium", "High"), lty=c(1,2,3), lwd=c(3,3,3), 
       col = c("green", "red", "blue"))
dev.off()

c_mtd = survfit(survie ~ m_td, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CHF_EXP_II_MASKER_TONE_DURATION_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mtd, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
     main = "", xlab = "Time (sec)", ylab = "Cumulative Hazard Function", las=1, 
     cex.axis=1.3, cex.lab=1.3, fun="cumhaz")
legend("topleft", c("Short", "Intermediate", "Long"), lty=c(1,2,3), 
       lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

c_tra = survfit(survie ~ t_ra, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CHF_EXP_II_TARGET_RATE_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_tra, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
     main = "", xlab = "Time (sec)", ylab = "Cumulative Hazard Function", las=1, 
     cex.axis=1.3, cex.lab=1.3, fun="cumhaz")
legend("topleft", c("5Hz","10Hz","20Hz"), lty=c(1,2,3), 
       lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

c_diss = survfit(survie ~ Diss, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CHF_EXP_II_MASKER_TARGET_SIMILARITY_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_diss, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3), 
     main = "", xlab = "Time (sec)", ylab = "Cumulative Hazard Function", las=1, 
     cex.axis=1.3, cex.lab=1.3, fun="cumhaz")
legend("topleft", c("High Similarity", "Medium Similarity", "Low Similarity"), 
       lty=c(1,2,3), lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

################################################################################
# CUMULATIVE DISTRIBUTION FUNCTION EN NOIR ET BLANC ET TRAIT POINTILLE
################################################################################
c_mppo = survfit(survie ~ m_ppo, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CDF_EXP_II_MASKER_SPECTRAL_DENSITY_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mppo, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Cumulative Distribution Function", las=1, cex.axis=1.3, 
     cex.lab=1.3, fun="F")
legend("topleft", c("Low", "Medium", "High"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_mtd = survfit(survie ~ m_td, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CDF_EXP_II_MASKER_TONE_DURATION_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mtd, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Cumulative Distribution Function", las=1, cex.axis=1.3, 
     cex.lab=1.3, fun="F")
legend("topleft", c("Short","Intermediate","Long"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_tra = survfit(survie ~ t_ra, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CDF_EXP_II_TARGET_RATE_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_tra, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab = "Cumulative Distribution Function", las=1, cex.axis=1.3, 
     cex.lab=1.3, fun="F")
legend("topleft", c("5Hz","10Hz","20Hz"), lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

c_diss = survfit(survie ~ Diss, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CDF_EXP_II_MASKER_TARGET_SIMILARITY_NB.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_diss, lty=c(1,2,3), lwd=c(3,3,3), main = "", xlab = "Time (sec)", 
     ylab="Cumulative Distribution Function", las=1, cex.axis=1.3, cex.lab=1.3, 
     fun="F")
legend("topleft", c("High Similarity", "Medium Similarity", "Low Similarity"), 
       lty=c(1,2,3), lwd=c(3,3,3))
dev.off()

################################################################################
# CUMULATIVE DISTRIBUTION FUNCTION EN COULEUR ET TRAIT POINTILLE
################################################################################
c_mppo = survfit(survie ~ m_ppo, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CDF_EXP_II_MASKER_SPECTRAL_DENSITY_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mppo, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3), 
     main = "", xlab = "Time (sec)", ylab = "Cumulative Distribution Function", 
     las=1, cex.axis=1.3, cex.lab=1.3, fun="F")
legend("topleft", c("Low", "Medium", "High"), lty=c(1,2,3), lwd=c(3,3,3), 
       col = c("green", "red", "blue"))
dev.off()

c_mtd = survfit(survie ~ m_td, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CDF_EXP_II_MASKER_TONE_DURATION_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_mtd, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
     main = "", xlab = "Time (sec)", ylab = "Cumulative Distribution Function", 
     las=1, cex.axis=1.3, cex.lab=1.3, fun="F")
legend("topleft", c("Short", "Intermediate", "Long"), lty=c(1,2,3), 
       lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

c_tra = survfit(survie ~ t_ra, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CDF_EXP_II_TARGET_RATE_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_tra, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3),  
     main = "", xlab = "Time (sec)", ylab = "Cumulative Distribution Function", 
     las=1, cex.axis=1.3, cex.lab=1.3, fun="F")
legend("topleft", c("5Hz","10Hz","20Hz"), lty=c(1,2,3), 
       lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

c_diss = survfit(survie ~ Diss, data = surv_df)
figure_to_save = "/home/link/Documents/thèse_onera/articles_alex/article_1_PA/v21/figures/figure_courbes_frailtymodel_all/CDF_EXP_II_MASKER_TARGET_SIMILARITY_COL.jpeg"
jpeg(file = figure_to_save, width=656, height=477, units = "px")
par(mar=c(4,4,0,0)+0.2)
plot(c_diss, col = c("green", "red", "blue"), lty=c(1,2,3), lwd=c(3,3,3), 
     main = "", xlab = "Time (sec)", ylab = "Cumulative Distribution Function", 
     las=1, cex.axis=1.3, cex.lab=1.3, fun="F")
legend("topleft", c("High Similarity", "Medium Similarity", "Low Similarity"), 
       lty=c(1,2,3), lwd=c(3,3,3), col = c("green", "red", "blue"))
dev.off()

################################################################################
