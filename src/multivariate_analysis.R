# Multivariate analysis of biochemical assayss

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

setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/")

# Here I would read in the three univariate files, but I already have them loaded. Remember for future runs.

PODWide <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/PeroxidaselDataWideFormat25Oct2017.csv")

ProteinWide <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/ProteinWideFormat25Oct2017.csv")

PPOWide <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/PolyphenolDataWideFormat25Oct2017.csv")


BiocAssaysWide <- merge(ProteinWide, PODWide)
BiocAssaysWide <- merge(BiocAssaysWide, PPOWide)

write.table(BiocAssaysWide, file = "AllBiocAssaysMANOVADataStyle14Nov2017.csv", sep = ",", row = F)

# Time for the big jump! 
BCmod1 <- lm(cbind(ProteinHr0, ProteinHr4, ProteinHr24, PODHr0, PODHr4, PODHr24,PPOHr0, PPOHr4, PPOHr24) ~ Range, data = BiocAssaysWide)

Assay <- factor(rep(c("Protein", "POD", "PPO"), each = 3))
MeasTime <- factor(rep(c("Hr0", "Hr4", "Hr24"), 3), levels = c("Hr0", "Hr4", "Hr24"))


BCmod2 <- lm(cbind(ProteinHr0, ProteinHr4, ProteinHr24, PODHr0, PODHr4, PODHr24,PPOHr0, PPOHr4, PPOHr24) ~ Site, data = BiocAssaysWide)

idata <- data.frame(Assay, MeasTime)

anova.BCmod1 <- Anova(BCmod1, idata= idata, idesign = ~ Assay * MeasTime, iterm = "Assay:MeasTime", type = 3)

anova.BCmod2 <- Anova(BCmod2, idata= idata, idesign = ~ Assay * MeasTime, type = 3)



pdf("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/results/exploratory/TempHE.pdf")
par(mfrow=c(3,3))
for(i in 1:9){
  for(j in 1:9){
    heplot(BCmod2, variables = c(i,j), cex = 0.5)
    
  }
}
dev.off()

heplot(BCmod1, idata= idata, idesign = ~ Assay * MeasTime, iterm = "Assay:MeasTime")
heplot(BCmod2, idata= idata, idesign = ~ Assay * MeasTime, iterm = "Assay:MeasTime")

CorData <- BiocAssaysWide[, c(6:14)]

CDmatrix1 <- cor(CorData)
CDmatrix2 <- rcorr(as.matrix(CorData))
CDmatrix3 <- rcorr.adjust(as.matrix(CorData))

corrplot(CDmatrix1, type = "upper")

corrplot(CDmatrix2$r, type = "upper", p.mat = CDmatrix2$P, sig.level = 0.05, insig = "blank")

corrplot(CDmatrix3$R$r, type = "upper", p.mat = apply(CDmatrix3$P, 2, as.numeric), sig.level = 0.05, insig = "blank")

PvalueAdj <- as.matrix(CDmatrix3$P)
PmatNum <- apply(PvalueAdj, 2, as.numeric)



#corrplot(CDmatrix$R$r, type = "upper", p.mat = CDmatrix$P.unadjust, sig.level = 0.05, insig = "blank", method = "number")


# ===== Making a figure of the model Nov 14, 2017 ======
# Need to use original data that has not been log transformed
ProteinUntra <- read.csv("/Users/chanj/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/ProteinData22Feb2017.csv")

PODUntra <- read.csv("/Users/chanj/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/WSU_PODfiles11Sept2017.csv")

PPOUntra <- read.csv("/Users/chanj/Documents/Friesen\ lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/WSU_PPOfiles11Sept2017.csv")

# Remove all extraneous columns
ProteinUnneeded <- names(ProteinUntra) %in% c("Well", "Weight", "Hour", "PCheck", "Plate", "fileName", "Absorbance")
ProteinUntra <- ProteinUntra[!ProteinUnneeded]

PODKeep <- names(PODUntra) %in% c("Sample", "Replicate", "Time", "AbsFreshWeight", "Range", "Genotype", "Site")
PODUntra <- PODUntra[PODKeep]

PPOKeep <- names(PPOUntra) %in% c("Sample", "Replicate", "Time", "AbsFreshWeight", "Range", "Genotype", "Site")
PPOUntra <- PPOUntra[PPOKeep]

# Rename response variable so that can bind the rows
names(ProteinUntra)[names(ProteinUntra) == "Protein"] <- "Response"
names(PODUntra)[names(PODUntra) == "AbsFreshWeight"] <- "Response"
names(PPOUntra)[names(PPOUntra) == "AbsFreshWeight"] <- "Response"

# Add column that identifies each assay
ProteinUntra$Assay <- rep("Protein", nrow(ProteinUntra))
PODUntra$Assay <- rep("POD", nrow(PODUntra))
PPOUntra$Assay <- rep("PPO", nrow(PPOUntra))

# Combine datasets
AllBiocUntrans <- rbind(ProteinUntra, PODUntra)
AllBiocUntrans <- rbind(AllBiocUntrans, PPOUntra)

# Make the fig!
ggplot(AllBiocUntrans, aes(x = Time, y = Response, shape = Range, colour = Range)) + stat_summary(fun.data = "mean_se") + facet_wrap(~Assay, scales = "free") + stat_summary(fun.y = mean, geom = "line", aes(group = Range))
ABU <- na.omit(AllBiocUntrans)
ggplot(ABU, aes(x = as.factor(Time), y = Response, shape = Range, colour = Range)) + stat_summary(fun.data = "mean_se") + facet_wrap(~Assay, scales = "free") + stat_summary(fun.y = mean, geom = "line", aes(group = Range))

# =============================

# ======= Centering MANOVA data Nov 14, 2017========
center_scale <- function(x) {
  scale(x, scale = FALSE)
}
BiocAssaysWideCenter <- center_scale(BiocAssaysWide[, 6:14])
BiocAssaysWide <- cbind(BiocAssaysWide, BiocAssaysWideCenter)
BAWC <- BiocAssaysWide[, -c(6:14)]

BCmod3 <- lm(cbind(ProteinHr0, ProteinHr4, ProteinHr24, PODHr0, PODHr4, PODHr24,PPOHr0, PPOHr4, PPOHr24) ~ Range, data = BAWC)

Assay <- factor(rep(c("Protein", "POD", "PPO"), each = 3))
MeasTime <- factor(rep(c("Hr0", "Hr4", "Hr24"), 3), levels = c("Hr0", "Hr4", "Hr24"))

BCmod4 <- lm(cbind(ProteinHr0, ProteinHr4, ProteinHr24, PODHr0, PODHr4, PODHr24,PPOHr0, PPOHr4, PPOHr24) ~ Site, data = BAWC)

idata <- data.frame(Assay, MeasTime)

anova.BCmod3 <- Anova(BCmod3, idata= idata, idesign = ~ Assay * MeasTime, iterm = "Assay:MeasTime", type = 3)

anova.BCmod4 <- Anova(BCmod4, idata= idata, idesign = ~ Assay * MeasTime, iterm = "Assay:MeasTime", type = 3)

# ====== Number of each population and by Range Nov 14, 2017 =====
BAWC %>% 
  group_by(Site) %>%
  summarise(no_rows = length(Site))

BAWC %>% count(Range, Site)
PopulationTable <- BAWC %>% count(Range, Site)
PopulationTable$NoGenotypes <- PopulationTable$n/3
PopTab <- dcast(PopulationTable, Site ~ Range, value.var = "NoGenotypes")

ggplot(PopulationTable, aes(Site, n/3, fill = Range, label = n/3)) + geom_bar(stat = "identity") + ylab("Number of Genotypes") + theme(axis.text.x = element_text(angle = 45)) + geom_text()


# ==== getting average values for each population for each assay
Pro0 <- ddply(ProteinUntra, c("Site", "Hour"),summarise, MeanPro0 = mean(Protein))
ProWide <- dcast(Pro0, Site ~ Hour)
PrefTab <- read.csv("Data/ProcessedData/ChoiceExpts/PrefTablePoints.csv")
PrefTab$ProteinHr0 <- ProWide[match(PrefTab$X, ProWide$Site), "0hr"]
PrefTab2 <- subset(PrefTab, !is.na(PrefTab$ProteinHr0))
cor(PrefTab2$SLMatchPoints, PrefTab2$ProteinHr0)
ggplot(PrefTab2, aes(x= ProteinHr0, y = VMatchPoints, colour = Range)) + geom_point() 
cor.test(PrefTab2$VMatchPoints, PrefTab2$ProteinHr0)

PrefTab2$Range <- ProteinWide[match(PrefTab2$X, ProteinWide$Site), "Range"]
mod <- lm(SLMatchPoints ~ Range, data = PrefTab2)
anova(mod)
ggplot(PrefTab2, aes(x= ProteinHr0, y = SLMatchPoints, colour = Range)) + geom_point() 
