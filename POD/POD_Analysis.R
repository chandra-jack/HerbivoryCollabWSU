# Feb 23, 2017
# Beginning analysis of protein quantification data

setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/POD")

library(platetools)
library(plater)
library(plyr)
library(dtplyr) # new combination of data.table and dplyr  ** Remember to install !!!!
library(lme4)
library(effects)
library(multcomp)
library(lmerTest)
# Prepping files: Treatment Plates
DataFiles <- read_plates(files = c("1POD_4hr_Plate3.csv", "2POD_4hr_Plate2.csv",  "4POD_24hr_Plate3.csv", "6POD_24hr_Plate2.csv", "7POD_24hr_Plate1.csv", "8POD_4hr_Plate1.csv", "9POD_0hr_Plate3.csv", "10POD_0hr_Plate2.csv", "11POD_0hr_Plate1.csv"), well_ids_column = "Well")

DataFiles<- rename(DataFiles, fileName = Plate)

DataFiles$Plate <- sapply(strsplit(DataFiles$fileName, "_"), "[",3)
DataFiles$Hour <- sapply(strsplit(DataFiles$fileName, "_"), "[",2)

filenames <- list.files(pattern = glob2rx("PODPlate*.csv"))
MapFiles <- lapply(filenames, function (x){
  dat <- read.csv(x)
  parts = strsplit(x, "_")[[1]]
  dat$Hour = as.factor(parts[2])
  dat$PCheck = as.factor(parts[3])
  return(dat)
})

MF <- bind_rows(MapFiles)
MF$Plate <- sapply(strsplit(MF$PCheck, ".csv"), "[",1)

PODdataFiles <- inner_join(MF, DataFiles, by = c("Well", "Plate", "Hour"))
PODdataFiles$MatchValue <- interaction(PODdataFiles$Sample, PODdataFiles$Hour)


# Prepping files: Control Plates

ControlFiles <- read_plates(files = c("3POD_4hr_Control.csv", "5POD_0hr_Control.csv", "12POD_24hr_Control.csv"), well_ids_column = "Well")

ControlFiles<- rename(ControlFiles, fileName = Plate)
ControlFiles$Plate <- sapply(strsplit(ControlFiles$fileName, "_"), "[",3)
ControlFiles$Hour <- sapply(strsplit(ControlFiles$fileName, "_"), "[",2)

MF2 <- subset(MF, Plate == "Plate4")
ControldataFiles <- inner_join(MF2, ControlFiles, by = c("Well", "Hour"))
ControldataFiles$MatchValue <- interaction(ControldataFiles$Sample, ControldataFiles$Hour)

# Adding Control data to treatment data
PODdataFiles$CtrlVal <- ControldataFiles[match(PODdataFiles$MatchValue, ControldataFiles$MatchValue), 10]

PODdataFiles$ABStrt_ABSctr <- PODdataFiles$values - PODdataFiles$CtrlVal
PODdataFiles$AbsFreshWeight <- PODdataFiles$ABStrt_ABSctr / PODdataFiles$Weight
PODdataFiles$Sample <- gsub("O", "0", PODdataFiles$Sample)
ProteinFiles <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/ProteinQuantification/ProteinData22Feb2017.csv")

PODdataFiles$Genotype <- ProteinFiles[match(PODdataFiles$Sample, ProteinFiles$Sample), 13]
PODdataFiles$Range <- ProteinFiles[match(PODdataFiles$Sample, ProteinFiles$Sample), 12]
PODdataFiles$Site <- ProteinFiles[match(PODdataFiles$Sample, ProteinFiles$Sample), 14]

PODdataFiles$Hour <- factor(PODdataFiles$Hour, levels = c("0hr", "4hr", "24hr"))

write.table(PODdataFiles, file = "WSU_PODfiles11Sept2017.csv", sep = ",", row = F)

ggplot(subset(PODdataFiles, Plate == "Plate1"), aes(Range, AbsFreshWeight)) + geom_boxplot() + facet_wrap( ~ Hour, scales = "free") + geom_hline(yintercept = 0, colour = "red") + labs(title ="Plate 1 POD data")
ggsave(file = "POD_plate1_23Feb2017.pdf")

ggplot(subset(PODdataFiles, Plate == "Plate2"), aes(Range, AbsFreshWeight)) +geom_boxplot() + facet_wrap( ~ Hour, scales = "free") + geom_hline(yintercept = 0, colour = "red") + labs(title ="Plate 2 POD data")
ggsave(file = "POD_plate2_23Feb2017.pdf")

ggplot(subset(PODdataFiles, Plate == "Plate3"), aes(Range, AbsFreshWeight)) +geom_boxplot() + facet_wrap( ~ Hour, scales = "free") + geom_hline(yintercept = 0, colour = "red") + labs(title ="Plate 3 POD data")
ggsave(file = "POD_plate3_23Feb2017.pdf")



SampleLists <- unique(PODdataFiles$Site)
plot_list = list()
for (i in seq_along(SampleLists)) {
  p = ggplot(subset(PODdataFiles, PODdataFiles$Site == SampleLists[i]), aes(x = Time, y = AbsFreshWeight, colour = Range)) +stat_summary(fun.data = "mean_se", size = 0.75) + stat_summary(geom="line", fun.y="mean")  + facet_wrap(~ Sample) + ggtitle(as.character(SampleLists[i]))
  plot_list[[i]] = p
}

pdf("PODBySite23Feb2017.pdf")
for (i in seq_along(SampleLists)) {
  print(plot_list[[i]])
}
dev.off()



ggplot(PODdataFiles, aes(Hour, AbsFreshWeight)) + stat_summary(fun.data = "mean_se", size = 0.5) + facet_wrap( ~ Range)  + labs(title = "POD data")
ggsave(file = "WSU_POD_28July2017.pdf")

pod_mod1 <- lm(AbsFreshWeight ~ Range * Hour * Site, data = PODdataFiles)
anova(pod_mod1)

# Cursory scan of the PODdataFiles set, shows some missing data that needs to be looked at more closely to fix.

PODsumm <- ddply(PODdataFiles, c("Range", "Site", "Hour"), summarise, Mean = mean(AbsFreshWeight, na.rm = T))

# make data wide
POD_wide<- dcast(PODsumm, Range + Site ~ Hour, value.var = "Mean")
POD_wide$Diff0_24 <- POD_wide$`24hr` - POD_wide$`0hr`

ggplot(PODsumm, aes(y = Site, x = Hour)) + geom_tile(aes(fill = Mean)) + facet_wrap(~ Range, scales = "free")
ggsave(file = "PODheat.pdf", width = 10, height = 8)
