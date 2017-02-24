# Feb 19, 2017
# Beginning analysis of protein quantification data
setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/ProteinQuantification")
library(platetools)
library(plater)
library(plyr)

# Prepping files
DataFiles <- read_plates(files = c("Plate1_0hr_17Feb2017.csv", "Plate2_0hr_17Feb2017.csv", "Plate3_0hr_17Feb2017.csv", "Plate1_4hr_17Feb2017.csv", "Plate2_4hr_17Feb2017.csv", "Plate3_4hr_17Feb2017.csv", "Plate1_24hr_17Feb2017.csv", "Plate2_24hr_17Feb2017.csv", "Plate3_24hr_17Feb2017.csv"), well_ids_column = "Well")

DataFiles<- rename(DataFiles, c("Plate" = "fileName"))

DataFiles$Plate <- sapply(strsplit(DataFiles$fileName, "_"), "[",1)
DataFiles$Hour <- sapply(strsplit(DataFiles$fileName, "_"), "[",2)

#parts = strsplit(DataFiles$Plate, "_")[[1]]

filenames <- list.files(pattern = glob2rx("ProteinQuant*.csv"))
MapFiles <- lapply(filenames, function (x){
  dat <- read.csv(x)
  parts = strsplit(x, "_")[[1]]
  dat$Hour = as.factor(parts[2])
  dat$PCheck = as.factor(parts[3])
  return(dat)
})

MF <- bind_rows(MapFiles)
MF$Plate <- sapply(strsplit(MF$PCheck, ".csv"), "[",1)

#MF$Plate <- rep(c("Plate1", "Plate2", "Plate3"), each = 96)

ProteinFiles <- inner_join(MF, DataFiles, by = c("Well", "Plate", "Hour"))
ProteinFiles<- rename(ProteinFiles, c("values" = "Absorbance"))

StCData <- read_plate(file = "StandardCurve_17Feb2017.csv", well_ids_column = "Well")
StdMap <- read.csv("ProteinStandardCurvePlateLayout.csv")
StandardCurve <- inner_join(StCData, StdMap, by = c("Well"))
StandardCurve <- rename(StandardCurve, c("values" = "x", "Sample" = "y"))
p <-ggplot(StandardCurve, aes(x, y)) + geom_point() + geom_smooth(method = "lm")
m <- lm(y ~ x, data = StandardCurve)
eq <- substitute(italic(y)== a + b %.% italic(x)*","~~italic(r)^2~"="~r2, list(a = format(coef(m)[1], digits = 4), b = format(coef(m)[2], digits = 4), r2 = format(summary(m)$r.squared, digits = 3)))
dftext <- data.frame(x = 0.5, y = 1500, eq = as.character(as.expression(eq)))
p + geom_text(aes(label = eq), data = dftext, parse = TRUE) + labs(y = "BSA Concentration (ug/mL", x = "Absorbance [562nm]", title = "Protein Standard Curve 17 Feb, 2017")
ggsave(file = "ProteinStandardCurve17Feb2017.pdf")

# Getting Protein Values
ProteinFiles <- within(ProteinFiles, Protein <- (-60.62 + 957.9 * Absorbance * 10)/ Weight)

ggplot(ProteinFiles, aes(x = Protein, colour = Range)) + geom_density() + facet_wrap(~ Hour, scales = "free_x")
ggplot(ProteinFiles, aes(Range, log(Protein), fill = Hour)) + geom_boxplot()
ProteinFiles$Hour <- factor(ProteinFiles$Hour, levels = c("0hr", "4hr", "24hr"))
ddply(ProteinFiles, c("Hour", "Range"), summarise, PMedian = median(Protein, na.rm = T), Pmean = mean(Protein, na.rm = T))

# To add which range the plant is from
InvadedVec <- c("W0076", "W0517", "W0130", "W0603", "W0146", "W0607", "W0153", "W0610", "W0155", "W0156", "W0167")

ProteinFiles$Range <- ifelse(ProteinFiles$Sample %in% InvadedVec, "Invasive", "Native")

# Trying to add Population data
MPmaster <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Copy of WSUV_MpolMlup_BulkingDesign_3.2016_masterfile -1.csv")
# Need to fix the mistaken Unique ID: WO to W0
ProteinFiles$Sample <- gsub("O", "0", ProteinFiles$Sample)
ProteinFiles$Genotype <- MPmaster[match(ProteinFiles$Sample, MPmaster$Unique.ID), 11]
ProteinFiles$Site <- MPmaster[match(ProteinFiles$Sample, MPmaster$Unique.ID), 1]

# Need to use as.character to prevent R from putting the level of the factor

which(is.na(ProteinFiles$Site))
ProteinFiles$Sample[19]      
ProteinFiles$Site <- ifelse(ProteinFiles$Genotype == "PI478440", "Bol", as.character(ProteinFiles$Site)) 

# Made a mistake in template files where I didn't change 0 to 24 for time; this fixes
ProteinFiles$Time <- ifelse(ProteinFiles$Hour == "24hr", 24, as.integer(ProteinFiles$Time))

ProteinFiles$Sample <- as.factor(ProteinFiles$Sample)

SampleLists <- unique(ProteinFiles$Site)
plot_list = list()
for (i in seq_along(SampleLists)) {
  p = ggplot(subset(ProteinFiles, ProteinFiles$Site == SampleLists[i]), aes(x = Time, y = Protein)) +stat_summary(fun.data = "mean_se", size = 0.75) + stat_summary(geom="line", fun.y="mean")  + facet_wrap(~ Sample) + ggtitle(as.character(SampleLists[i]))
  plot_list[[i]] = p
}

 
#ggsave(p, filename = paste("Site_",SampleLists[i],"22Feb2017", ".pdf", sep = ""))

pdf("ProteinQuantificationBySite23Feb2017.pdf")
for (i in seq_along(SampleLists)) {
  print(plot_list[[i]])
}
dev.off()

#write.table(ProteinFiles, file = "ProteinData22Feb2017.csv", sep = ",", row = )
ProteinFiles <- read.csv("ProteinData22Feb2017.csv")
