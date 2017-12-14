 # Dec 7, 2016

# Beginning data analysis of WSU collab project

# First starting with Preference 1 Data worksheet


library(ggplot2)
library(data.table)
library(scales)
library(plyr)
library(reshape2)
library(lme4)
library(effects)
library(multcomp)
library(lmerTest)
library(piecewiseSEM)
library(car)
library(gridExtra)
library(cowplot)

setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/")

ConstData <- read.csv("Cynthia_Data/ConstPref_28Apr2017.csv")

# There seems to be a mistake in the data labeling for V004, so will remove for now.

ConstData <- subset(ConstData, AssayID != "V004")

# Start with some exploratory graphics
# Win or Lose- where a value of 1 means it was preferred
ggplot(ConstData, aes(x = Range, y = WL_own))+ stat_summary(fun.data = "mean_se", size = 0.75) + geom_hline(aes(yintercept = 0.5), colour = "red") + theme_classic() + theme(axis.ticks.x = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12), strip.background = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12)) + labs(y = "Winning %" ) + scale_y_continuous(labels = percent) + facet_wrap(~InsectSpp)

ggplot(ConstData, aes(x = Range, y = WL_total))+ stat_summary(fun.data = "mean_se", size = 0.75) + geom_hline(aes(yintercept = 0.5), colour = "red") + theme_classic() + theme(axis.ticks.x = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12), strip.background = element_blank(), strip.text.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12)) + labs(y = "Winning %" ) + scale_y_continuous(labels = percent) + facet_wrap(~InsectSpp)

ggplot(VelvetbeanData, aes(x = Site, y = WL_own, colour = Range))+ stat_summary(fun.data = "mean_se", size = 0.75) + geom_hline(aes(yintercept = 0.5), colour = "red") + theme_classic() + theme(axis.ticks.x = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12), strip.background = element_blank(), strip.text.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12)) + labs(y = "Winning %" ) + scale_y_continuous(labels = percent) + facet_wrap(~Range, scales = "free_x")

ggplot(ConstData, aes(x = Site, y = WL_total, colour = Range))+ stat_summary(fun.data = "mean_se", size = 0.75) + geom_hline(aes(yintercept = 0.5), colour = "red") + theme_classic() + theme(axis.ticks.x = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12), strip.background = element_blank(), strip.text.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12)) + labs(y = "Winning %" ) + scale_y_continuous(labels = percent) + facet_grid(InsectSpp~Range, scales = "free_x")

ConstData2 <- ConstData
ConstData2$WL_own <- as.factor(ConstData2$WL_own)
mod <- glm(WL_own ~ Range/Site + InsectSpp, family = binomial(link = "logit"), data = ConstData2)
anova(mod, test = "Chisq")

mod <- glmer(WL_own ~ Range +InsectSpp + (1|Site), family = binomial(link = "logit") , data = ConstData2)
anova(mod, test = "Chisq")
summary(mod)
anova(mod)

ConstData2$WL_total <- as.factor(ConstData2$WL_total)
mod <- glm(WL_total ~ Range/Site + InsectSpp, family = binomial(link = "logit"), data = ConstData2)
anova(mod, test = "Chisq")

mod <- glmer(WL_total ~ Range + InsectSpp+ (1|Site), family = binomial(link = "logit") , data = ConstData2)
anova(mod, test = "Chisq")
summary(mod)
anova(mod)


ggplot(ConstData, aes(x = Range, y = Ratio_Consumed))+ stat_summary(fun.data = "mean_se", size = 0.75) + geom_hline(aes(yintercept = 1), colour = "red") + theme_classic() + theme(axis.ticks.x = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12), strip.background = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12)) + facet_wrap(~InsectSpp, scales = "free_x")

ggplot(VelvetbeanData, aes(x = Range, y = Ratio_Consumed))+ stat_summary(fun.data = "mean_se", size = 0.75) + geom_hline(aes(yintercept = 1), colour = "red") + theme_classic() + theme(axis.ticks.x = element_blank(), legend.position = "bottom", legend.text = element_text(size = 12), strip.background = element_blank(), strip.text.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12)) + labs(y = "Winning %" ) 
