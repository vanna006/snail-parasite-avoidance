########################################
# R code for Physa movement experiment #
########################################

setwd('E:/A_in_arev/Physa_move')
shake <- read.csv('shake.csv')

shake1 <- subset(shake, day == '1d')
shake3 <- subset(shake, day == '3d')
shake5 <- subset(shake, day == '5d')
shake3.5 <- rbind(shake3,shake5)

#### let's look at shaking first
#day 1
hist(shake1$shakes, breaks=15)

library(glmmTMB)
fit_hzip <- glmmTMB(shakes ~ group + (1|arena), data=shake1, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)

#test alternate models for shakes
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom2, family=list(family="truncated_nbinom1",link="log"))

library(bbmle)
AICtab(fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)

#the best one
summary(fit_zinbinom2)

#day3
fit_hzip <- glmmTMB(shakes ~ group + (1|arena), data=shake3, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)

#test alternate models for shakes
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom1, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)
#the best one
summary(fit_hnbinom2)

#day 3 and 5 (day 5 did not work alone)
hist(shake3.5$shakes, breaks=15)
fit_hzip <- glmmTMB(shakes ~ group + (1|arena) + (1|snail), data=shake3.5, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)

#test alternate models for shakes
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom2, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)
#the best one
summary(fit_hnbinom2)


#### look at surfacing behavior
#day1
hist(shake1$surfacing, breaks = 15)
hist(shake3.5$surfacing, breaks = 15)

library(lme4)
mm1 <- glmer.nb(surfacing ~ group + (1|arena), data = shake1)
summary(mm1)

#day3.5
mm2 <- glmer.nb(surfacing ~ group + (1|snail), data = shake3.5)
summary(mm2)


#### look at time eating
hist(shake1$time.eating, breaks=15)

mm3 <- glmer(time.eating ~ group + length + (1|arena), data = shake1, family = 'poisson')
summary(mm3)

plot(time.eating ~ length, col = group, pch = 16, data = shake1)

#day3.5
mm4 <- glmer(time.eating ~ group + length + (1|arena) +(1|snail), data = shake3.5, family = 'poisson')
summary(mm4)

plot(time.eating ~ length, col = group, pch = 16, data = shake3.5)

#as a function of infection intensity
#mm5 <- glmer.nb(time.eating ~ metas + (1|arena) + (1|snail), data = shake3.5)
#summary(mm5)


#### look at mean velocity
hist(shake1$mean.vel, breaks = 15)
shapiro.test(shake1$mean.vel)

mm5 <- lmer(mean.vel ~ group + (1|arena), data = shake1)
summary(mm5)
library(car)
Betas  <- fixef(mm5)                  #Get the betas
SE     <-  sqrt(diag(vcov(mm5)))      #Get the SEs
pval   <- 2*pnorm(-abs(Betas  / SE)) #Z distribution
Output <- cbind(Betas,SE, pval)
print(Output, digits = 3) # P values for fixed effect
plot(mm5)
library(lattice)
qqmath(mm5)
Anova(mm5, type="III")

hist(shake3.5$mean.vel, breaks = 15)
shapiro.test(shake3.5$mean.vel)
#looks gamma to me
mm6 <- glmer(mean.vel ~ group + (1|arena) + (1|snail), data = shake3.5, family = Gamma(link = "inverse"))
summary(mm6)
qqmath(mm6)

#### look at variance in velocity
hist(shake1$vel.var, breaks = 15)
shapiro.test(shake1$vel.var)
mm7 <- lmer(vel.var ~ group + (1|arena), data = shake1)
summary(mm7)
qqmath(mm7)
Anova(mm7, type="III")

hist(shake3.5$vel.var, breaks = 15)
shapiro.test(shake3.5$vel.var^0.4044118)
bestNormalize(shake3.5$vel.var)

mm8 <- lmer((vel.var^0.4044118) ~ group + (1|arena) + (1|snail), data = shake3.5)
summary(mm8)
qqmath(mm8)
Anova(mm8, type="III")


#### look at variance in turn direction
hist(shake1$turn.var, breaks = 15)
shapiro.test(shake1$turn.var)
mm9 <- lmer(turn.var ~ group + (1|arena), data = shake1)
summary(mm9)
qqmath(mm9)
Anova(mm9, type="III")

hist(shake3.5$turn.var, breaks = 15)
shapiro.test(shake3.5$turn.var)
mm10 <- lmer(turn.var ~ group + (1|arena), data = shake3.5)
summary(mm10)
qqmath(mm10)
Anova(mm10, type="III")


#### let's see if anything influences meta intensity
hist(shake1$metas)
shapiro.test(shake1$metas)

mm11 <- lmer(metas ~ shakes + (1|arena), data = shake1)
summary(mm11)
qqmath(mm11)
Anova(mm11, type="III")

mm12 <- lmer(metas ~ bouts + (1|arena), data = shake1)
summary(mm12)
qqmath(mm12)
Anova(mm12, type="III")

mm13 <- lmer(metas ~ shake.bout + (1|arena), data = shake1)
summary(mm13)
qqmath(mm13)
Anova(mm13, type="III")

mm14 <- lmer(metas ~ surfacing + (1|arena), data = shake1)
summary(mm14)
qqmath(mm14)
Anova(mm14, type="III")

mm15 <- lmer(metas ~ mean.vel + (1|arena), data = shake1)
summary(mm15)
qqmath(mm15)
Anova(mm15, type="III")

mm16 <- lmer(metas ~ vel.var + (1|arena), data = shake1)
summary(mm16)
qqmath(mm16)
Anova(mm16, type="III")

mm17 <- lmer(metas ~ turn.var + (1|arena), data = shake1)
summary(mm17)
qqmath(mm17)
Anova(mm17, type="III")

mm18 <- lmer(metas ~ time.eating + (1|arena), data = shake1)
summary(mm18)
qqmath(mm18)
Anova(mm18, type="III")

plot(metas ~ time.eating, data = shake1)

library(ggplot2)
shakeplot <- ggplot(shake1, aes(x = group, y = shakes, fill = group)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Shakes per trial", hjust=10) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))

saltlabels <- c('Unexposed', 'Exposed')

a <- shakeplot + geom_jitter(size = 3, color = 'grey25', width = 0.15, height=0) + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 1, position=position_dodge(1), color='white') +
  scale_fill_manual(values=c("white","grey75")) + theme_classic() +
  theme(text = element_text(size=12), legend.position='none') + scale_x_discrete(labels= saltlabels)


surfplot <- ggplot(shake1, aes(x = group, y = surfacing, fill = group)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Surfacing events per trial", hjust=10) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))


b <- surfplot + geom_jitter(size = 3, color = 'grey25', width = 0.15, height = 0) + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 1, position=position_dodge(1), color = 'white') +
  scale_fill_manual(values=c("white","grey75")) + theme_classic() +
  theme(text = element_text(size=12), legend.position='none') + scale_x_discrete(labels= saltlabels)


eatplot <- ggplot(shake1, aes(x = group, y = time.eating, fill = group)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Time spent foraging (seconds per trial)", hjust=10) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))


c <- eatplot + geom_jitter(size = 3, color = 'grey25', width = 0.15, height = 0) + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 1, position=position_dodge(1), color='white') +
  scale_fill_manual(values=c("white","grey75")) + theme_classic() +
  theme(text = element_text(size=12), legend.position='none') + scale_x_discrete(labels= saltlabels)


speedplot <- ggplot(shake1, aes(x = group, y = mean.vel, fill = group)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Mean velocity (mm/s)", hjust=10) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))


d <- speedplot + geom_jitter(size = 3, color = 'grey25', width = 0.15, height = 0) + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 1, position=position_dodge(1), color='white') +
  scale_fill_manual(values=c("white","grey75")) + theme_classic() +
  theme(text = element_text(size=12), legend.position='none') + scale_x_discrete(labels= saltlabels)


speedvar <- ggplot(shake1, aes(x = group, y = vel.var, fill = group)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Variance in velocity (mm/s)", hjust=10) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))


speedvar + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("white","grey75")) + theme_classic() +
  theme(text = element_text(size=10), legend.position='none') + scale_x_discrete(labels= saltlabels)

library(ggpubr)
ggarrange(a, b, c, d , 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

