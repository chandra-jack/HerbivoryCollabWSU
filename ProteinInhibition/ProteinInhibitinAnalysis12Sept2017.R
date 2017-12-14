# WSU samples Proteinase Inhibition

library(ggplot2)
library(reshape2)
library(platetools)
library(plater)
library(plyr)
library(dplyr)
library(data.table)

setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/ProteinInhibition")

# Combine plate data files
ListPIData <- list.files(pattern = glob2rx("*hr.csv"))
PIdataFiles <- read_plates(ListPIData, well_ids_column = "Well")

PIdataFiles <- rename(PIdataFiles, fileName = Plate)
PIdataFiles$Plate <- sapply(strsplit(PIdataFiles$fileName, "_"), "[",2)
PIdataFiles$Hour <- sapply(strsplit(PIdataFiles$fileName, "_"), "[",3)

# Combine plate description files
ListPIDesc <- list.files(pattern = glob2rx("*A.csv"))
MapFiles <- lapply(ListPIDesc, function (x){
  dat <- read.csv(x)
  parts = strsplit(x, "_")[[1]]
  dat$Hour = as.factor(parts[3])
  dat$Plate = as.factor(parts[4])
  return(dat)
})

MF <- bind_rows(MapFiles)

# Join data and description files
PIfiles <- inner_join(MF, PIdataFiles, by = c("Well", "Plate", "Hour"))
# Forgot to remove the empty data wells
PIfiles <- subset(PIfiles, Sample != "")
PIfiles$Time <- ifelse(PIfiles$Hour == "24hr", 24, PIfiles$Time)
# Calculations for PI
# Subtract TRIS sample from Trypsin sample
PIfiles$SubAbs <- PIfiles$Abs_A - PIfiles$Abs_B

# Separate Samples from Controls
PIfiles$Sample <- as.factor(PIfiles$Sample)
PIfiles$Hour <- as.factor(PIfiles$Hour)

PI_Cont <- subset(PIfiles, Sample == "Control") 
PI_C2 <- ddply(PI_Cont, "Hour", summarise, meanAbs = mean(SubAbs), seAbs = sqrt(var(SubAbs)/length(SubAbs)))

# Jump to bottom of script (line: ) to do statistical analysis

PID <- ddply(PIfiles, c("Sample", "Time", "Weight", "Hour"), summarise, meanSubAbs = mean(SubAbs), seSubAbs = sqrt(var(SubAbs)/length(SubAbs)))


PI_Donly <- subset(PID, Sample != "Control")

PItry <- inner_join(PI_Donly, PI_C2, by = ("Hour"))

PItry <- within(PItry, {
  PI <- (1- (meanSubAbs/meanAbs)) * 100
  PIse1 <- abs((sqrt((seSubAbs/meanSubAbs)^2 + (seAbs/meanAbs)^2)/100) * (meanSubAbs/meanAbs))
  PIse2 <- abs(PIse1 * PI)
  PI_g <- PI/ Weight 
  PIse <- abs(PI_g * PIse2)
})

InvadedVec <- c("W0076", "W0517", "W0130", "W0603", "W0146", "W0607", "W0153", "W0610", "W0155", "W0156", "W0167")

PItry$Range <- ifelse(PItry$Sample %in% InvadedVec, "Invasive", "Native")
PItry$Sample <- gsub("O", "0", PItry$Sample)


PItry$Genotype <- ProteinFiles[match(PItry$Sample, ProteinFiles$Sample), "Genotype"]
PItry$Site <- ProteinFiles[match(PItry$Sample, ProteinFiles$Sample), "Site"]
write.table(PItry, file = "WSU_ProteinaseInhibition12Sept2017_absSE.csv",  sep = ",", row = F)

ProteinFiles$Sample <- as.factor(ProteinFiles$Sample)

SampleLists <- unique(PItry$Site)
plot_list = list()
for (i in seq_along(SampleLists)) {
  p = ggplot(subset(PItry, PItry$Site == SampleLists[i]), aes(x = Time, y = PI_g, colour = Range)) +stat_summary(fun.data = "mean_se", size = 0.75) + stat_summary(geom="line", fun.y="mean")  + facet_wrap(~ Sample) + ggtitle(as.character(SampleLists[i]))
  plot_list[[i]] = p
}

pdf("ProteinaseInhibitionBySite12Sept2017.pdf")
for (i in seq_along(SampleLists)) {
  print(plot_list[[i]])
}
dev.off()


# Needed to do some stats: 
PI_FilesDonly <- subset(PIfiles, Sample != "Control")

PItry2 <- inner_join(PI_FilesDonly, PI_C2, by = ("Hour"))
PItry2$PI <- (1- (PItry2$SubAbs/PItry2$meanAbs)) * 100

# The next section is for the methods paper 10 Oct 2017
ProIf <- subset(PItry2, subset = Sample %in% c("W0419" ,"W0420", "W0077", "W0607", "W0079", "W0076", "W0517", "W0603", "W0146", "W0421"))
ProIf <- subset(ProIf, Time != 4)

mod1 <- lm(PI ~ Time * Sample, data = ProIf)
anova(mod1)
ddply(ProIf,"Time", summarise, MEAN = mean(PI), SEval = sqrt(var(PI)/length(PI)))
