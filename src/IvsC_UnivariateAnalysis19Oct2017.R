# Univariate Analyses for the Biochemical Assays of the IvsC project
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

# Data clean up
 # I considered combining data from all three assays and then cleaning but I'm just going to do it individually. I will combine when I get to multivariate analysis.

# ===== Oct 19, 2017 ====
setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/UnivariateAnalyses/")

# The protein, pod, and ppo data files all have similar structure

# ====Protein Quantification ====
ProteinFile <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/ProteinData22Feb2017.csv")

# Removed unnecessary columns
ProRemove <- names(ProteinFile) %in% c("PCheck", "Plate", "fileName", "Absorbance")
ProteinFile <- ProteinFile[!ProRemove]

# Removed USDA genotypes
MPmaster <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Copy of WSUV_MpolMlup_BulkingDesign_3.2016_masterfile -1.csv")

MP_USDA <- MPmaster$Unique.ID[MPmaster$Pod.Produced != "USDA"]

ProteinFile <- ProteinFile[ProteinFile$Sample %in% MP_USDA,]

ProteinFile$Hour <- factor(ProteinFile$Hour, levels = c("0hr", "4hr", "24hr"))

with(ProteinFile, interaction.plot(Hour, Range, Protein))
boxplot(Protein~Hour, data = ProteinFile)
ggplot(ProteinFile, aes(Hour, Protein, fill = Range)) + geom_boxplot() +ylim(0, 1*10^6)# + facet_wrap(~ Site, scales = "free")

# ==== List of outliers ====
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

dat <- ProteinFile %>% tibble::rownames_to_column(var="outlier") %>% group_by(Range, Hour) %>% mutate(is_outlier=ifelse(is_outlier(Protein), Protein, as.numeric(NA)))
dat$outlier[which(is.na(dat$is_outlier))] <- as.numeric(NA)
ggplot(dat, aes(y=Protein, x=Hour, fill= Range)) + geom_boxplot() + geom_text(aes(label=outlier),na.rm=TRUE,nudge_y=0.05)

ProtOutlier <- subset(dat, !is.na(outlier))

P1 <-ggplot(ProteinFile, aes(Hour, Protein, fill = Range)) + geom_boxplot()
P2 <- ggplot(ProteinFile, aes(Hour, Protein, fill = Range)) + geom_boxplot() +ylim(0, 1*10^6)
ProBoxplot <-ggdraw() + draw_plot(P1 + theme(legend.position = "none"), 0, 0, 1, 1) + draw_plot(P2, 0.5, 0.6, 0.5, 0.5)
save_plot("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/ProteinBoxplots19Oct2017.pdf", ProBoxplot, base_width = 11, base_height = 8)
ggplot(aes(Protein), data = ProteinFile) + geom_density() + facet_grid(Range ~ Hour)
ggsave("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/ProteinDensityFacetOct2017.pdf", width = 11, height = 8)

# ==== Transforming to normalize =====
ProteinFile$ProSqrt <- sqrt(ProteinFile$Protein)
qqnorm(ProteinFile$ProSqrt)
qqline(ProteinFile$ProSqrt) # Did not change much
densityPlot(ProteinFile$ProSqrt)


ProteinFile$ProCube <- ProteinFile$Protein ^ 0.33
qqnorm(ProteinFile$ProCube)
qqline(ProteinFile$ProCube)
densityPlot(ProteinFile$ProCube)

ProteinFile$ProLog <- log(ProteinFile$Protein)
qqnorm(ProteinFile$ProLog)
qqline(ProteinFile$ProLog)
densityPlot(ProteinFile$ProLog)


qqnorm(ProteinFile$Protein)
qqline(ProteinFile$Protein)
densityPlot(ProteinFile$Protein)

pdf("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/ProteinHistTrans19Oct2017.pdf")
par(mfrow=c(2,2))
densityPlot(ProteinFile$Protein)
densityPlot(ProteinFile$ProSqrt)
densityPlot(ProteinFile$ProCube)
densityPlot(ProteinFile$ProLog)
dev.off()

# will use log transformation


# Because this is a repeated measures (Time) model, can use a MANOVA for this
ProteinWide <- dcast(ProteinFile, Sample + Site + Replicate + Genotype + Range ~ Hour, value.var = "ProLog")
ProteinWide <- na.exclude(ProteinWide)
colnames(ProteinWide)[colnames(ProteinWide) =="0hr"] <- "ProteinHr0"
colnames(ProteinWide)[colnames(ProteinWide) =="4hr"] <- "ProteinHr4"
colnames(ProteinWide)[colnames(ProteinWide) =="24hr"] <- "ProteinHr24"
write.table(ProteinWide, file = "~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/ProteinWideFormat25Oct2017.csv", sep = ",", row = F)
ProMod1 <- lm(cbind(ProteinHr0, ProteinHr4, ProteinHr24) ~ Range, data = ProteinWide)
plot(density(resid(ProMod1)))
qqnorm(resid(ProMod1))
qqline(resid(ProMod1))

ProResd <- resid(ProMod1)

blockD <- ordered(c("Hr0", "Hr4", "Hr24"), levels = c("Hr0", "Hr4", "Hr24"))
contrasts(blockD) <- matrix(c(-1,1,0,0,-1,1), ncol = 2)
idata <-data.frame(blockD)

rm(blockD2)
proMod1a <-Anova(ProMod1, idata = idata, idesign = ~ blockD)
summary(Anova(ProMod1, idata = idata, idesign = ~ blockD, type = "III"))

heplot(ProMod1, idata = idata, idesign = ~ blockD, iterm =  "blockD", cex = 0.5)
heplot(ProMod1, variables = c(2,3))
pairs(ProMod1,cex = 0.5)

# exploratory graphs
PRdf <- data.frame(ProResd)
ggplot(PRdf) + stat_qq(aes(sample = ProteinHr0,colour = "blue")) + stat_qq(aes(sample = ProteinHr4,colour = "red")) + stat_qq(aes(sample = ProteinHr24,colour = "green")) + abline(y=1)+ geom_abline(intercept = 0, slope = 0.6)+ scale_colour_manual(name = 'legend',values =c('blue'='blue','red'='red', 'green' = 'green'), labels = c('Hour 0','Hour 4', 'Hour 24'))

f <- fitted(ProMod1)
r <- rstandard(ProMod1)

ggsave(file = "~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/ProteinResidual24Oct2017.pdf")

write.table(PODWide, file = "~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/ProteinDataWideFormat25Oct2017.csv", sep = ",", row = F)
# ====POD assay ==== Oct 21, 2017
PODfile <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/WSU_PODfiles11Sept2017.csv")

# Removed unnecessary columns
PODRemove <- names(PODfile) %in% c("PCheck", "Plate", "fileName", "values", "MatchValue", "CtrlVal","ABStrt_ABSctr")
PODfile <- PODfile[!PODRemove]

PODfile <- PODfile[PODfile$Sample %in% MP_USDA,]
PODfile$Hour <- factor(PODfile$Hour, levels = c("0hr", "4hr", "24hr"))
Pd1<- ggplot(PODfile, aes(Hour, AbsFreshWeight, fill = Range)) + geom_boxplot() 
Pd2 <- ggplot(PODfile, aes(Hour, AbsFreshWeight, fill = Range)) + geom_boxplot() +ylim(0, 1*10^2)
ProBoxplot <-
  ggdraw() + draw_plot(Pd1 + theme(legend.position = "none"), 0, 0, 1, 1) + draw_plot(Pd2, 0.1, 0.5, 0.5, 0.5)
save_plot("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/PeroxidaseBoxplots23Oct2017.pdf", ProBoxplot, base_width = 11, base_height = 8)
# seems to be a lot more variance in the bigger (Native) group

# checking normality

# standardized fresh weight values go negative
min(PODfile$AbsFreshWeight) # -3.75
# will add 4 so all will be positive for transformations
PODfile$AbsFreshWeightN <- PODfile$AbsFreshWeight + 4

qqnorm(PODfile$AbsFreshWeightN)
qqline(PODfile$AbsFreshWeightN)
densityPlot(PODfile$AbsFreshWeightN)

PODfile$PodSqrt <- sqrt(PODfile$AbsFreshWeightN)
qqnorm(PODfile$PodSqrt)
qqline(PODfile$PodSqrt) 
densityPlot(PODfile$PodSqrt)


PODfile$PodCube <- PODfile$AbsFreshWeightN ^ 0.33
qqnorm(PODfile$PodCube)
qqline(PODfile$PodCube)
densityPlot(PODfile$PodCube)

PODfile$PodLog <- log(PODfile$AbsFreshWeightN)
qqnorm(PODfile$PodLog)
qqline(PODfile$PodLog)
densityPlot(PODfile$PodLog)

F1 <-ggplot(PODfile, aes(PodLog, fill = Range)) + geom_density()  + facet_grid(Range~ Hour, scales = "free")
F2 <-ggplot(PODfile, aes(PodLog, fill = Range)) + geom_density()  + facet_grid(Range~ Hour)
PODdist <-plot_grid(F1, F2)

save_plot("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/PeroxidaseDist23Oct2017.pdf", PODdist, base_width = 11, base_height = 8)

# Because this is a repeated measures (Time) model, can use a MANOVA for this
PODWide <- dcast(PODfile, Sample + Site + Replicate + Genotype + Range ~ Hour, value.var = "PodLog")
PODWide <- na.exclude(PODWide)
colnames(PODWide)[colnames(PODWide) =="0hr"] <- "PODHr0"

colnames(PODWide)[colnames(PODWide) =="4hr"] <- "PODHr4"
colnames(PODWide)[colnames(PODWide) =="24hr"] <- "PODHr24"

PodMod1 <- lm(cbind(PODHr0, PODHr4, PODHr24) ~ Range, data = PODWide)
plot(density(resid(PodMod1)))
qqnorm(resid(PodMod1))
qqline(resid(PodMod1))

PodResd <- resid(PodMod1)

blockPOD <- ordered(c("PODHr0", "PODHr4", "PODHr24"), levels = c("PODHr0", "PODHr4", "PODHr24"))
contrasts(blockPOD) <- matrix(c(-1,1,0,0,-1,1), ncol = 2)
idata <- data.frame(blockPOD)

Anova(PodMod1, idata = idata, idesign = ~ blockPOD)
summary(Anova(PodMod1, idata = idata, idesign = ~ blockD), multivariate = F)

heplot(PodMod1, idata = idata, idesign = ~blockD, iterm = "blockD")

PODdf <- data.frame(PodResd)
ggplot(PODdf) + stat_qq(aes(sample = PODHr0,colour = "blue")) + stat_qq(aes(sample = PODHr4,colour = "red")) + stat_qq(aes(sample = PODHr24,colour = "green")) + scale_colour_manual(name = 'legend',values =c('blue'='blue','red'='red', 'green' = 'green'), labels = c('Hour 0','Hour 4', 'Hour 24'))

f2 <- fitted(PodMod1)
r2 <- rstandard(PodMod1)
plot(f2, r2, col = as.numeric(col(f2)), pch = 19)
legend("center", legend = paste0("response ", 1:ncol(r2)), pch = 19,col = 1:ncol(r2), text.col = 1:ncol(r2))

write.table(PODWide, file = "~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/PeroxidaselDataWideFormat25Oct2017.csv", sep = ",", row = F)



# ===== PPO =====

PPOFile <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/WSU_PPOfiles11Sept2017.csv")
PPOFile <- subset(PPOFile, Sample !="")

# Removed unnecessary columns
ProRemove <- names(PPOFile) %in% c("PCheck", "Plate", "fileName", "values", "MatchValue", "CtrlVal","ABStrt_ABSctr")
PPOFile <- PPOFile[!ProRemove]

# Removed USDA genotypes
MPmaster <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Copy of WSUV_MpolMlup_BulkingDesign_3.2016_masterfile -1.csv")

MP_USDA <- MPmaster$Unique.ID[MPmaster$Pod.Produced != "USDA"]

PPOFile <- PPOFile[PPOFile$Sample %in% MP_USDA,]
PPOFile$Hour <- factor(PPOFile$Hour, levels = c("0hr", "4hr", "24hr"))

ggplot(PPOFile, aes(Hour, AbsFreshWeight, fill = Range)) + geom_boxplot() 

ggsave("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/PolyphenolOBoxplots23Oct2017.pdf",width = 11, height = 8)

min(PPOFile$AbsFreshWeight) # -0.74
# will add 1 so all will be positive for transformations
PPOFile$AbsFreshWeightN <- PPOFile$AbsFreshWeight + 1

qqnorm(PPOFile$AbsFreshWeightN)
qqline(PPOFile$AbsFreshWeightN)
densityPlot(PPOFile$AbsFreshWeightN)

PPOFile$PpoSqrt <- sqrt(PPOFile$AbsFreshWeightN)
qqnorm(PPOFile$PpoSqrt)
qqline(PPOFile$PpoSqrt) 
densityPlot(PPOFile$PpoSqrt)


PPOFile$PpoCube <- PPOFile$AbsFreshWeightN ^ 0.33
qqnorm(PPOFile$PpoCube)
qqline(PPOFile$PpoCube)
densityPlot(PPOFile$PpoCube)

PPOFile$PpoLog <- log(PPOFile$AbsFreshWeightN)
qqnorm(PPOFile$PpoLog)
qqline(PPOFile$PpoLog)
densityPlot(PPOFile$PpoLog)

ggplot(PPOFile, aes(PpoLog, fill = Range)) + geom_density()  + facet_grid(Range~ Hour)
ggsave("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/PolyphenolODensity23Oct2017.pdf",width = 11, height = 8)

# Because this is a repeated measures (Time) model, can use a MANOVA for this
PPOWide <- dcast(PPOFile, Sample + Site + Replicate + Genotype + Range ~ Hour, value.var = "PpoLog")
PPOWide <- na.exclude(PPOWide)
colnames(PPOWide)[colnames(PPOWide) =="0hr"] <- "PPOHr0"

colnames(PPOWide)[colnames(PPOWide) =="4hr"] <- "PPOHr4"
colnames(PPOWide)[colnames(PPOWide) =="24hr"] <- "PPOHr24"

PPOMod1 <- lm(cbind(PPOHr0, PPOHr4, PPOHr24) ~ Range, data = PPOWide)
plot(density(resid(PPOMod1)))
qqnorm(resid(PPOMod1))
qqline(resid(PPOMod1))

blockPPO <- ordered(c("PPOHr0", "PPOHr4", "PPOHr24"), levels = c("PPOHr0", "PPOHr4", "PPOHr24"))
contrasts(blockPPO) <- matrix(c(-1,1,0,0,-1,1), ncol = 2)
idata <- data.frame(blockPPO)

Anova(PPOMod1, idata = idata, idesign = ~ blockPPO)
summary(Anova(PPOMod1, idata = idata, idesign = ~ blockPPO))

heplot(PPOMod1, cex = 0.6)
heplot(PPOMod1, idata = idata, idesign = ~blockPPO, iterm = "blockPPO")




PpoResd <- resid(PPOMod1)
PPOdf <- data.frame(PpoResd)
ggplot(PPOdf) + stat_qq(aes(sample = PPOHr0,colour = "blue")) + stat_qq(aes(sample = PPOHr4,colour = "red")) + stat_qq(aes(sample = PPOHr24,colour = "green")) + scale_colour_manual(name = 'legend',values =c('blue'='blue','red'='red', 'green' = 'green'), labels = c('Hour 0','Hour 4', 'Hour 24'))

f2 <- fitted(PPOMod1)
r2 <- rstandard(PPOMod1)
plot(f2, r2, col = as.numeric(col(f2)), pch = 19)
legend("topleft", legend = paste0("response ", 1:ncol(r2)), pch = 19,col = 1:ncol(r2), text.col = 1:ncol(r2))

write.table(PPOWide, file = "~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/PolyphenolDataWideFormat25Oct2017.csv", sep = ",", row = F)
#write.table(ProteinFiles, file = "ProteinData22Feb2017.csv", sep = ",", row = F)