# Protein Quantification Assay for WSU samples
# Plate Layout
# Feb 15, 2017


library(ggplot2)
library(platetools)
library(plater)
library(dplyr)
library(ggplot2bdc)
library(cowplot)

setwd("ProteinQuantification/")
Tissue1_2 <- read.csv("Tissue Collection for Chemical Analysis- Selected Plants - Plate #1 & #2.csv")
Tissue1_2 <- subset(Tissue1_2, Well.Number != "B2")
TFdf <- Tissue1_2[rep(row.names(Tissue1_2), each = 3),]
write.table(TFdf, file = "Plates1and2_replicateIDs.csv", sep = ",", row = F)

# ----------------------------- Code to test if cell is within a vector; 
# Plan to use this in conjunction with dataset that list whether invasive or not
n = c(2, 3, 5) 
s = c("aa", "bb", "cc") 
b = c(TRUE, FALSE, TRUE) 
df = data.frame(n, s, b) 

MyVec <- c("aa", "dd", "ee", "cc")
df$VecCK <- ifelse(df$s %in% MyVec, "YES", "NO")
# ----------------------------- Make plate layouts 

PQplatemap1 <- read.csv("ProteinQuantification/ProteinQuantPlateLayout_0hr_1of3.csv")
# Next line is needed to get the perfect looking plate layout
PQplatemap1 <- subset(PQplatemap1, Sample != "")
PQplatemap1 <- mutate(PQplatemap1, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

 pl1 <- ggplot(data=PQplatemap1, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 0 hour 1 of 3")
 
 PQplatemap2 <- read.csv("ProteinQuantification/ProteinQuantPlateLayout_0hr_2of3.csv")
 PQplatemap2 <- subset(PQplatemap2, Sample != "")
 PQplatemap2 <- mutate(PQplatemap2, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))
 
 pl2 <- ggplot(data=PQplatemap2, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 0 hour 2 of 3")
 
 PQplatemap3 <- read.csv("ProteinQuantPlateLayout_0hr_3of3.csv")
 PQplatemap3 <- subset(PQplatemap3, Sample != "")
 PQplatemap3 <- mutate(PQplatemap3, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))
 
 pl3 <- ggplot(data=PQplatemap3, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 0 hour 3 of 3")
 
 StdCurve <- read.csv("ProteinStandardCurvePlateLayout.csv")
 StdCurve <- subset(StdCurve, Sample != "")
 StdCurve <- mutate(StdCurve, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))
 
 stc <- ggplot(data=StdCurve, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(aes(colour = SampleType),size=8) +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 3) + labs(title = "Protein Quantification Standard Curve Layout")
 
 pgrid <- plot_grid(stc, pl1, nrow = 2)
 save_plot("ProteinPlateLabel_1of2Feb16_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)
 
 pgrid <- plot_grid(pl2, pl3, nrow = 2)
 save_plot("ProteinPlateLabel_2of2Feb16_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)
 
 
 
 # ---------------------------- 4 hour plate layout
 Tissue3_4 <- read.csv("Tissue Collection for Chemical Analysis- Selected Plants - Plate #3 & #4.csv")
 Tissue3_4 <- subset(Tissue3_4, Well.Number != "B2")
 TFdf2 <- Tissue3_4[rep(row.names(Tissue3_4), each = 3),]
 write.table(TFdf2, file = "Plates3and4_replicateIDs.csv", sep = ",", row = F)
 
 PQplatemap1_4 <- read.csv("ProteinQuantPlateLayout_4hr_1of3.csv")
 # Next line is needed to get the perfect looking plate layout
 PQplatemap1_4 <- subset(PQplatemap1_4, Sample != "")
 PQplatemap1_4 <- mutate(PQplatemap1_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))
 
pl1.4 <-ggplot(data=PQplatemap1_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 4 hour 1 of 3")


PQplatemap2_4 <- read.csv("ProteinQuantPlateLayout_4hr_2of3.csv")
# Next line is needed to get the perfect looking plate layout
PQplatemap2_4 <- subset(PQplatemap2_4, Sample != "")
PQplatemap2_4 <- mutate(PQplatemap2_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl2.4 <-ggplot(data=PQplatemap2_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 4 hour 2 of 3")


PQplatemap3_4 <- read.csv("ProteinQuantPlateLayout_4hr_3of3.csv")
PQplatemap3_4 <- subset(PQplatemap3_4, Sample != "")
PQplatemap3_4 <- mutate(PQplatemap3_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl3.4 <-ggplot(data=PQplatemap3_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") + coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) +  scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter()+ geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 4 hour 3 of 3")

pgrid <- plot_grid(pl1.4, pl2.4, nrow = 2)
save_plot("ProteinPlateLabel_4hr_1of2Feb17_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)

pgrid <- plot_grid(pl3.4,nrow = 2)
save_plot("ProteinPlateLabel_4hr_2of2Feb17_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)


# --------------------------- 24 hour plate layout
Tissue5_6 <- read.csv("Tissue Collection for Chemical Analysis- Selected Plants - Plate #5 & #6.csv")
Tissue5_6 <- subset(Tissue5_6, Well.Number != "B2")
TFdf3 <- Tissue5_6[rep(row.names(Tissue5_6), each = 3),]
write.table(TFdf3, file = "Plates5and6_replicateIDs.csv", sep = ",", row = F)

PQplatemap1_24 <- read.csv("ProteinQuantPlateLayout_24hr_1of3.csv")
# Next line is needed to get the perfect looking plate layout
PQplatemap1_24 <- subset(PQplatemap1_24, Sample != "")
PQplatemap1_24 <- mutate(PQplatemap1_24, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl1.4 <-ggplot(data=PQplatemap1_24, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 24 hour 1 of 3")


PQplatemap2_24 <- read.csv("ProteinQuantPlateLayout_24hr_2of3.csv")
# Next line is needed to get the perfect looking plate layout
PQplatemap2_24 <- subset(PQplatemap2_24, Sample != "")
PQplatemap2_24 <- mutate(PQplatemap2_24, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl2.4 <-ggplot(data=PQplatemap2_24, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 24 hour 2 of 3")


PQplatemap3_24 <- read.csv("ProteinQuantPlateLayout_24hr_3of3.csv")
PQplatemap3_24 <- subset(PQplatemap3_24, Sample != "")
PQplatemap3_24 <- mutate(PQplatemap3_24, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl3.4 <-ggplot(data=PQplatemap3_24, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") + coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) +  scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter()+ geom_text(aes(label = Sample), size = 2.1) + labs(title = "Protein Quantification Map 24 hour 3 of 3")

pgrid <- plot_grid(pl1.4, pl2.4, nrow = 2)
save_plot("ProteinPlateLabel_24hr_1of2Feb17_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)

pgrid <- plot_grid(pl3.4,nrow = 2)
save_plot("ProteinPlateLabel_24hr_2of2Feb17_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)

