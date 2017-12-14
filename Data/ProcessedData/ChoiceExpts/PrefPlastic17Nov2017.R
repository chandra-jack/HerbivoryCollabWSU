library(ggplot2)
library(data.table)
library(scales)
library(plyr)
library(dplyr)
library(reshape2)
library(lme4)
library(effects)
library(multcomp)
library(lmerTest)
library(piecewiseSEM)
library(car)
library(gridExtra)
library(cowplot)
library(Rmisc)
library(heplots)
library(Hmisc)
library(RcmdrMisc)
library(corrplot)


setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/ChoiceExpts/")

PrefData <- read.csv("Corrected_FINAL_PreferenceData.csv")
MPmaster <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/WSUV_BulkingData_Final_ZCL.csv")

# Preparing the dataset
PrefData$Genotype_Plant1 <- MPmaster[match(PrefData$PlantID_1, MPmaster$Unique.ID), "MSU.Genotype"]
PrefData$Population_Plant1 <- MPmaster[match(PrefData$PlantID_1, MPmaster$Unique.ID), "Site.ID"]

PrefData$Genotype_Plant2 <- MPmaster[match(PrefData$PlantID_2, MPmaster$Unique.ID), "MSU.Genotype"]
PrefData$Population_Plant2 <- MPmaster[match(PrefData$PlantID_2, MPmaster$Unique.ID), "Site.ID"]

#write.table(PrefData, file = "ModifiedPreferenceData16Nov2017.csv", sep = ",", row = F)

PrefData$TissuePref <- PrefData$AmtPlant1Consumed - PrefData$AmtPlant2Consumed

PrefData$TissuePrefScaled <- (PrefData$AmtPlant1Consumed - PrefData$AmtPlant2Consumed)/(PrefData$AmtPlant1Consumed + PrefData$AmtPlant2Consumed)
# Running a glm binomial model to check for preference for invasive or native population
mod1 <- glm(as.factor(Inv_1_Nat_0) ~ Herbivore*HerbWt_Final , family = binomial(link = "logit"), data = PrefData)
anova(mod1, test = "Chisq")
#summary(mod1)

ggplot(PrefData, aes(CatStage, Inv_1_Nat_0, colour = Herbivore)) + stat_summary(fun.data = "mean_se") + geom_hline(yintercept = 0.5, colour = "black")

# Separating herbivore data and then running one sample wilcox tests

PrefLoop <- subset(PrefData, Herbivore == "Soybean Looper")
PrefV <- subset(PrefData, Herbivore != "Soybean Looper")

wilcox.test(PrefLoop$Inv_1_Nat_0, mu = 0.5)
wilcox.test(PrefV$Inv_1_Nat_0, mu = 0.5)
kruskal.test(PrefLoop$Inv_1_Nat_0 ~ PrefLoop$CatStage)
# Comparison based on amount of tissue consumed, not just which was consumed more
ggplot(PrefData, aes(CatStage, TissuePref, colour = Herbivore)) + stat_summary(fun.data = "mean_se") + geom_hline(yintercept = 0, colour = "black") 

mod3 <- lm(TissuePref ~ Herbivore*HerbWt_Final, data = PrefData)
anova(mod3)
M3.res <- rstandard(mod3)
qqnorm(M3.res)
qqline(M3.res)

ggplot(PrefData, aes(CatStage, TissuePrefScaled, colour = Herbivore)) + stat_summary(fun.data = "mean_se") + geom_hline(yintercept = 0, colour = "black")

mod4 <- lm(TissuePrefScaled ~ Herbivore*HerbWt_Final, data = PrefData)
anova(mod4)
# Also ran t test
wilcox.test(PrefLoop$TissuePrefScaled, mu = 0)
wilcox.test(PrefLoop$TissuePref, mu = 0)
kruskal.test(PrefLoop$TissuePrefScaled ~ PrefLoop$CatStage)

wilcox.test(PrefV$TissuePrefScaled, mu = 0)
kruskal.test(PrefV$TissuePref ~ PrefV$CatStage)

#=====Plastic Data=====
PlasticData <- read.csv("Corrected_FINAL_PlasticityData.csv")

# Preparing the dataset
PlasticData$Genotype_Plant1 <- MPmaster[match(PlasticData$PlantID_1, MPmaster$Unique.ID), "MSU.Genotype"]
PlasticData$Population_Plant1 <- MPmaster[match(PlasticData$PlantID_1, MPmaster$Unique.ID), "Site.ID"]

PlasticData$Genotype_Plant2 <- MPmaster[match(PlasticData$PlantID_2, MPmaster$Unique.ID), "MSU.Genotype"]
PlasticData$Population_Plant2 <- MPmaster[match(PlasticData$PlantID_2, MPmaster$Unique.ID), "Site.ID"]

write.table(PlasticData, file = "ModifiedPlasticity7Dec2017cj.csv", sep = ",", row = F)
PlasticData$TissuePref <- PlasticData$AmtPlant1Consumed - PlasticData$AmtPlant2Consumed
mod1 <- glm(as.factor(Con_1_Ind_0) ~ Herbivore*HerbWt_Final + Range, family = binomial(link = "logit"), data = PlasticData)
anova(mod1, test = "Chisq")
#summary(mod1)

ggplot(PlasticData, aes(CatStage, Con_1_Ind_0, colour = Herbivore)) + stat_summary(fun.data = "mean_se") + geom_hline(yintercept = 0.5, colour = "black") + facet_wrap(~ Range)


PlasticLoop <- subset(PlasticData, Herbivore == "Soybean Looper")
PlasticV <- subset(PlasticData, Herbivore != "Soybean Looper")

wilcox.test(PlasticLoop$Con_1_Ind_0 ~ PlasticLoop$Range)
wilcox.test(PlasticV$Con_1_Ind_0 ~ PlasticV$Range)

# tissue amount 
mod4 <- lm(TissuePref ~ Herbivore*HerbWt_Final + Range, data = PlasticData)
anova(mod4)

#mod3 <- lm(TissuePref ~ Herbivore*CatStage + Population_Plant1:Population_Plant2, data = PlasticData)
#mod5 <- lm(TissuePref ~ Herbivore*CatStage + Population_Plant1 + Population_Plant2, data = PlasticData)
#AIC(mod3, mod4, mod5)

ggplot(PlasticData, aes(CatStage, TissuePref, colour = Herbivore)) + stat_summary(fun.data = "mean_se") + geom_hline(yintercept = 0, colour = "black") + facet_wrap(~ Range)

# Looper
wilcox.test(PlasticLoop$TissuePref ~ PlasticLoop$Range)
#Velvet bean
wilcox.test(PlasticV$TissuePref ~ PlasticV$Range)

# ==============================
# 7 Dec 2017  looking at cat weight and cat stage
wtMod <- lm(HerbWt_Final ~ Herbivore * CatStage, data = PlasticData)
anova(wtMod)


# going to try and match up bioc data
BiocSub2 <- ddply(BiocSub, c("Sample", "Assay", "Time"), summarise, meanResponse = mean(Response))
PopTab <- dcast(PopulationTable, Site ~ Range, value.var = "NoGenotypes")

CastBioc <- dcast(BiocSub2, Sample  ~ Assay +Time, value.var = "meanResponse")

PrefData$Population_Plant2 <- MPmaster[match(PrefData$PlantID_2, MPmaster$Unique.ID), "Site.ID"]

BiocPref <- PrefData
BiocPref$Protein0 <- CastBioc[match(BiocPref$PlantID_1, CastBioc$Sample), "Protein_0"]

BiocPref$Protein0 <- CastBioc[match(BiocPref$PlantID_2, CastBioc$Sample), 9]



# ======== Number of genotypes used =====
Genos <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/AllGenotypes.csv")
sapply(Genos, function(x) length(unique(x)))

Genos$Range <- MPmaster[match(Genos$Genotypes, MPmaster$Unique.ID), "Range"]

Genos$Site <- MPmaster[match(Genos$Genotypes, MPmaster$Unique.ID), "Origin"]

ddply(Genos, .(Range), 'nrow')
write.table(Genos, "Supplement.csv", sep = ",", row.names = FALSE)
