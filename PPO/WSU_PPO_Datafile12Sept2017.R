# WSU herbivory project- 11 Sept 2017: Polyphenoloxidase

library(platetools)
library(plater)
library(plyr)
library(dplyr)
library(ggplot2)
library(readxl) # allows reading in excel files


setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/PPO")

# POD and PPO use same plate layout template
# use the same files (copied over) and change name to PPO prefix
f <- list.files(pattern = glob2rx("POD*.csv"))
fout <- gsub("POD", "PPO", f)
for (i in 1:length(f)){
  file.rename(f[i], fout[i])
}

# Renaming the datafiles from .xls to .csv
#- First create a new directory to keep the .xls files for safekeeping
dir.create("ExcelDataFiles")
f_fin <- list.files(pattern = glob2rx("PPO_*.xlsx"))
file.copy(f_fin, "ExcelDataFiles")

# To convert files from xlsx to csv:
lapply(f_fin, function(f) {
  df = read_excel(f, sheet = 1)
  write.csv(df, gsub("xlsx", "csv", f), row.names = F)
})

# Used regexp to pull out just the data files separate from the controls
Dlist <- list.files(pattern = glob2rx("*hr_pl*.csv")) # data
Clist <- list.files(pattern = glob2rx("*hr_co*.csv")) # controls

# Read in all the plates
PPODataFiles <- read_plates(Dlist, well_ids_column = "Well")

# Rename Plate to filename for next manipulation
PPODataFiles<- rename(PPODataFiles, fileName = Plate)

PPODataFiles$Plate <- sapply(strsplit(PPODataFiles$fileName, "_"), "[",3)
PPODataFiles$Hour <- sapply(strsplit(PPODataFiles$fileName, "_"), "[",2)

filenames <- list.files(pattern = glob2rx("PPOPlate*.csv"))
MapFiles <- lapply(filenames, function (x){
  dat <- read.csv(x)
  parts = strsplit(x, "_")[[1]]
  dat$Hour = as.factor(parts[2])
  dat$PCheck = as.factor(parts[3])
  return(dat)
})

MF <- bind_rows(MapFiles)
MF$Plate <- sapply(strsplit(MF$PCheck, ".csv"), "[",1)

library(Hmisc) # plate is not capitalized in PPODataFiles and needs to be changed
PPODataFiles$Plate <- capitalize(PPODataFiles$Plate)
PPODataF <- inner_join(MF, PPODataFiles, by = c("Well", "Plate", "Hour"))
PPODataF$MatchValue <- interaction(PPODataF$Sample, PPODataF$Hour)

# Prepping Control Files
ControlFiles <- read_plates(Clist, well_ids_column = "Well")

ControlFiles<- rename(ControlFiles, fileName = Plate)
ControlFiles$Plate <- sapply(strsplit(ControlFiles$fileName, "_"), "[",3)
ControlFiles$Hour <- sapply(strsplit(ControlFiles$fileName, "_"), "[",2)

MF2 <- subset(MF, Plate == "Plate4")
ControldataFiles <- inner_join(MF2, ControlFiles, by = c("Well", "Hour"))
ControldataFiles$MatchValue <- interaction(ControldataFiles$Sample, ControldataFiles$Hour)

# Adding Control data to treatment data
PPODataF$CtrlVal <- ControldataFiles[match(PPODataF$MatchValue, ControldataFiles$MatchValue), 10]

PPODataF$ABStrt_ABSctr <- PPODataF$values - PPODataF$CtrlVal
PPODataF$AbsFreshWeight <- PPODataF$ABStrt_ABSctr / PPODataF$Weight
PPODataF$Sample <- gsub("O", "0", PPODataF$Sample)
ProteinFiles <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/ProteinQuantification/ProteinData22Feb2017.csv")

# Combining data and description files
PPODataF$Genotype <- ProteinFiles[match(PPODataF$Sample, ProteinFiles$Sample), 13]
PPODataF$Range <- ProteinFiles[match(PPODataF$Sample, ProteinFiles$Sample), 12]
PPODataF$Site <- ProteinFiles[match(PPODataF$Sample, ProteinFiles$Sample), 14]

PPODataF$Hour <- factor(PPODataF$Hour, levels = c("0hr", "4hr", "24hr"))

write.table(PPODataF, file = "WSU_PPOfiles11Sept2017.csv", sep = ",", row = F)

SampleLists <- unique(PPODataF$Site)
plot_list = list()
for (i in seq_along(SampleLists)) {
  p = ggplot(subset(PPODataF, PPODataF$Site == SampleLists[i]), aes(x = Time, y = Protein, colour = Range)) +stat_summary(fun.data = "mean_se", size = 0.75) + stat_summary(geom="line", fun.y="mean")  + facet_wrap(~ Sample) + ggtitle(as.character(SampleLists[i]))
  plot_list[[i]] = p
}