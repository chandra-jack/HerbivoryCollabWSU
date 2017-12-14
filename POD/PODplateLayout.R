# POD  Assay for WSU samples
# Plate Layout
# Feb 23, 2017


library(ggplot2)
library(platetools)
library(plater)
library(dplyr)
library(ggplot2bdc)
library(cowplot)


setwd("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/POD/")

PODplatemap1 <- read.csv("PODPlateLayout_0hr_Plate1.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap1 <- subset(PODplatemap1, Sample != "")
PODplatemap1 <- mutate(PODplatemap1, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl1.0 <-ggplot(data=PODplatemap1, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 0 hour 1 of 3")


PODplatemap2 <- read.csv("PODPlateLayout_0hr_Plate2.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap2 <- subset(PODplatemap2, Sample != "")
PODplatemap2 <- mutate(PODplatemap2, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl2.0 <-ggplot(data=PODplatemap2, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 0 hour 2 of 3")

PODplatemap3 <- read.csv("PODPlateLayout_0hr_Plate3.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap3 <- subset(PODplatemap3, Sample != "")
PODplatemap3 <- mutate(PODplatemap3, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl3.0 <-ggplot(data=PODplatemap3, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 0 hour 3 of 3")

PODplatemap4 <- read.csv("PODPlateLayout_0hr_Plate4.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap4 <- subset(PODplatemap4, Sample != "")
PODplatemap4 <- mutate(PODplatemap4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl4.0 <-ggplot(data=PODplatemap4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 0 hour CONTROL")


pgrid <- plot_grid(pl1.0, pl2.0, nrow = 2)
save_plot("PODPlateLabel_0hr_1of2Feb23_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)

pgrid <- plot_grid(pl3.0, pl4.0, nrow = 2)
save_plot("PODPlateLabel_0hr_2of2Feb23_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)




#------
PODplatemap1_4 <- read.csv("PODPlateLayout_4hr_Plate1.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap1_4 <- subset(PODplatemap1_4, Sample != "")
PODplatemap1_4 <- mutate(PODplatemap1_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl1.4 <-ggplot(data=PODplatemap1_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 4 hour 1 of 3")


PODplatemap2_4 <- read.csv("PODPlateLayout_4hr_Plate2.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap2_4 <- subset(PODplatemap2_4, Sample != "")
PODplatemap2_4 <- mutate(PODplatemap2_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl2.4 <-ggplot(data=PODplatemap2_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 4 hour 2 of 3")

PODplatemap3_4 <- read.csv("PODPlateLayout_4hr_Plate3.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap3_4 <- subset(PODplatemap3_4, Sample != "")
PODplatemap3_4 <- mutate(PODplatemap3_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl3.4 <-ggplot(data=PODplatemap3_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 4 hour 3 of 3")

PODplatemap4_4 <- read.csv("PODPlateLayout_4hr_Plate4.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap4_4 <- subset(PODplatemap4_4, Sample != "")
PODplatemap4_4 <- mutate(PODplatemap4_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl4.4 <-ggplot(data=PODplatemap4_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 4 hour CONTROL")


pgrid <- plot_grid(pl1.4, pl2.4, nrow = 2)
save_plot("PODPlateLabel_4hr_1of2Feb23_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)

pgrid <- plot_grid(pl3.4, pl4.4, nrow = 2)
save_plot("PODPlateLabel_4hr_2of2Feb23_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)



#------
PODplatemap1_4 <- read.csv("PODPlateLayout_4hr_Plate1.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap1_4 <- subset(PODplatemap1_4, Sample != "")
PODplatemap1_4 <- mutate(PODplatemap1_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl1.4 <-ggplot(data=PODplatemap1_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 4 hour 1 of 3")


PODplatemap2_4 <- read.csv("PODPlateLayout_4hr_Plate2.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap2_4 <- subset(PODplatemap2_4, Sample != "")
PODplatemap2_4 <- mutate(PODplatemap2_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl2.4 <-ggplot(data=PODplatemap2_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 4 hour 2 of 3")

PODplatemap3_4 <- read.csv("PODPlateLayout_4hr_Plate3.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap3_4 <- subset(PODplatemap3_4, Sample != "")
PODplatemap3_4 <- mutate(PODplatemap3_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl3.4 <-ggplot(data=PODplatemap3_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 4 hour 3 of 3")

PODplatemap4_4 <- read.csv("PODPlateLayout_4hr_Plate4.csv")
# Next line is needed to get the perfect looking plate layout
PODplatemap4_4 <- subset(PODplatemap4_4, Sample != "")
PODplatemap4_4 <- mutate(PODplatemap4_4, Row=as.numeric(match(toupper(substr(Well, 1, 1)), LETTERS)), Column=as.numeric(substr(Well, 2, 5)))

pl4.4 <-ggplot(data=PODplatemap4_4, aes(x=Column, y=Row)) + geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2), color="grey90", fill="white", shape=21, size=6) + geom_point(size=8, colour = "lightgreen") +  coord_fixed(ratio=(13/12)/(9/8), xlim = c(0.5, 12.5), ylim=c(0.5, 8.5)) + scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) + scale_x_continuous(breaks=seq(1, 12)) + theme_bdc_microtiter() + geom_text(aes(label = Sample), size = 2.1) + labs(title = "POD Map 4 hour CONTROL")


pgrid <- plot_grid(pl1.4, pl2.4, nrow = 2)
save_plot("PODPlateLabel_4hr_1of2Feb23_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)

pgrid <- plot_grid(pl3.4, pl4.4, nrow = 2)
save_plot("PODPlateLabel_4hr_2of2Feb23_2017.pdf", pgrid, nrow = 2, base_aspect_ratio = 1.5)
