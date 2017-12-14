

# Attempting to resample Bioc assay values
BiocAssaysWide <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/AllBiocAssaysMANOVADataStyle14Nov2017.csv")

BV_Nat <- subset(BiocAssaysWide, Range == "Native")
BV_Inv <- subset(BiocAssaysWide, Range == "Invasive")

# There are only 33 invasive -> because of replication, 11 pops
# Therefore we will need only 11 pops of native each iteration

# First method: randomly grab invasive populations. calculation the variances and determine if they differ from the variance of the full population

N = 1000
TestMat <- matrix(rep(0, 9*N), ncol = 9)
for (i in 1:N){
  pops <- sample(unique(BV_Nat$Genotype), 11) # Sample populations randomly 
  BV_Nat2 <- BV_Nat[BV_Nat$Genotype %in% pops, ] # Narrow dataset
  TestMat[i,] <- apply(BV_Nat2[,6:ncol(BV_Nat2)], 2, function(x) var(x, na.rm=TRUE))
}


VecVar <- apply(BV_Nat[,6:ncol(BV_Nat)], 2, function(x) var(x, na.rm=TRUE))
ColnNames <- colnames(BV_Nat)
ColnNames <- ColnNames[-(1:5)]
NewMat <- matrix(rep(0, 27), ncol = 3)
for (i in 1:9) {
  
  NewMat[i,1] <- ColnNames[i]
  NewMat[i,2] <- t.test(TestMat[,i], mu = VecVar[i])$p.value
  NewMat[i,3] <- mean(abs(TestMat[,i] -mean(TestMat[,i])) > VecVar[i])
}
NewMat

# Method 2: Resample native population, combine with invasive, run levene's test
N = 1000
TestMat <- matrix(rep(0, 9*N), ncol = 9)
for (i in 1:N){
  pops <- sample(unique(BV_Nat$Genotype), 11, replace = T) # Sample populations randomly 
  BV_Nat2 <- BV_Nat[BV_Nat$Genotype %in% pops, ] # Narrow dataset
  FullBV <- rbind(BV_Nat2, BV_Inv)
  test <- leveneTests(FullBV[,6:14], FullBV$Range)
  FStat <- test[,3]
  for (j in 1:9){
    TestMat[i, j] <- FStat[j]
  }
}

ColnNames <- colnames(BV_Nat)
ColnNames <- ColnNames[-(1:5)]
testObs <- leveneTests(BiocAssaysWide[,6:14], BiocAssaysWide$Range)
FStatObs <- testObs[,3]

FMat <- matrix(rep(0, 18), ncol = 2)
for (i in 1:9) {
  FMat[i,1] <- ColnNames[i]
  FMat[i,2] <- mean(abs(TestMat[,i] -mean(TestMat[,i])) > FStatObs[i])
}
FMat

PODWide <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/PeroxidaselDataWideFormat25Oct2017.csv")
PPOWide <- read.csv("~/Documents/Friesen lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/PolyphenolDataWideFormat25Oct2017.csv")
ProWide <- read.csv("/Users/chanj/Documents/Friesen\ lab/MedicagoHerbPopulation/HerbivoryCollabWSU/Data/ProcessedData/ProteinWideFormat25Oct2017.csv")

block <- ordered(c("Hr0", "Hr4", "Hr24"), levels = c("Hr0", "Hr4", "Hr24"))
contrasts(block) <- matrix(c(-1,1,0,0,-1,1), ncol = 2)
idata <- data.frame(block)

PODWide_Nat <- subset(PODWide, Range == "Native")
PODWide_Inv <- subset(PODWide, Range == "Invasive")

N = 1000
PodMat <- matrix(rep(0, 6*N), ncol = 6)

for (i in 1:N) {
  outtests <- car:::print.Anova.mlm
  body(outtests)[[16]] <- quote(invisible(tests))
  body(outtests)[[15]] <- NULL
  pops <- sample(unique(PODWide_Nat$Genotype), 11, replace = T) 
  PODWide_Nat2 <- PODWide_Nat[PODWide_Nat$Genotype %in% pops, ] # Narrow dataset
  FullPODWide <- rbind(PODWide_Nat2, PODWide_Inv)
  PodMod1 <- lm(cbind(PODHr0, PODHr4, PODHr24) ~ Range, data = FullPODWide)
  tab <- lapply(c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),  function(i)  outtests(Anova(PodMod1, idata = idata, idesign = ~ block, test.statistic=i)))
  tab <- do.call(rbind, tab)
  Fval <- tab[,3]
  Pval <- tab[,6]
  for (j in 1:3){
    X = j + 1
    Y = j + 3
    PodMat[i,j] <- Fval[X]
    PodMat[i,Y] <- Pval[X]
  }
}
  

PPOWide_Nat <- subset(PPOWide, Range == "Native")
PPOWide_Inv <- subset(PPOWide, Range == "Invasive")
N = 1000
PpoMat <- matrix(rep(0, 6*N), ncol = 6)

for (i in 1:N) {
  outtests <- car:::print.Anova.mlm
  body(outtests)[[16]] <- quote(invisible(tests))
  body(outtests)[[15]] <- NULL
  pops <- sample(unique(PPOWide_Nat$Genotype), 11, replace = T) 
  PPOWide_Nat2 <- PPOWide_Nat[PPOWide_Nat$Genotype %in% pops, ] # Narrow dataset
  FullPPOWide <- rbind(PPOWide_Nat2, PPOWide_Inv)
  PpoMod1 <- lm(cbind(PPOHr0, PPOHr4, PPOHr24) ~ Range, data = FullPPOWide)
  tab <- lapply(c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),  function(i)  outtests(Anova(PpoMod1, idata = idata, idesign = ~ block, test.statistic=i)))
  tab <- do.call(rbind, tab)
  Fval <- tab[,3]
  Pval <- tab[,6]
  for (j in 1:3){
    X = j + 1
    Y = j + 3
    PpoMat[i,j] <- Fval[X]
    PpoMat[i,Y] <- Pval[X]
  }
}


ProWide_Nat <- subset(ProWide, Range == "Native")
ProWide_Inv <- subset(ProWide, Range == "Invasive")
N = 1000
ProMat <- matrix(rep(0, 6*N), ncol = 6)

for (i in 1:N) {
  outtests <- car:::print.Anova.mlm
  body(outtests)[[16]] <- quote(invisible(tests))
  body(outtests)[[15]] <- NULL
  pops <- sample(unique(PODWide_Nat$Genotype), 11, replace = T) 
  ProWide_Nat2 <- ProWide_Nat[ProWide_Nat$Genotype %in% pops, ] # Narrow dataset
  FullProWide <- rbind(ProWide_Nat2, ProWide_Inv)
  ProMod1 <- lm(cbind(ProteinHr0, ProteinHr4, ProteinHr24) ~ Range, data = FullProWide)
  tab <- lapply(c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),  function(i)  outtests(Anova(ProMod1, idata = idata, idesign = ~ block, test.statistic=i)))
  tab <- do.call(rbind, tab)
  Fval <- tab[,3]
  Pval <- tab[,6]
  for (j in 1:3){
    X = j + 1
    Y = j + 3
    ProMat[i,j] <- Fval[X]
    ProMat[i,Y] <- Pval[X]
  }
}

 
# Run actual values from dataset

#pod
outtests <- car:::print.Anova.mlm
body(outtests)[[16]] <- quote(invisible(tests))
body(outtests)[[15]] <- NULL
PodMod1 <- lm(cbind(PODHr0, PODHr4, PODHr24) ~ Range, data = PODWide)
tab <- lapply(c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),  function(i)  outtests(Anova(PodMod1, idata = idata, idesign = ~ block, test.statistic=i)))
tab <- do.call(rbind, tab)

Fval <- tab[,3]
Pval <- tab[,6]
PODtrue <- data.frame(Fval[2:4], Pval[2:4])
rownames(PODtrue) <- c("Range", "block", "Range:block")
names(PODtrue)[names(PODtrue) == 'Fval.2.4.'] <- 'test.statistic'
names(PODtrue)[names(PODtrue) == 'Pval.2.4.'] <- 'p.value'

for (i in 1:3) {
  PODtrue$BootF[i] <-mean(abs(PodMat[,i] - mean(PodMat[,i])) > PODtrue$test.statistic[i])
}

# PPO
outtests <- car:::print.Anova.mlm
body(outtests)[[16]] <- quote(invisible(tests))
body(outtests)[[15]] <- NULL
PpoMod1 <- lm(cbind(PPOHr0, PPOHr4, PPOHr24) ~ Range, data = PPOWide)
tab <- lapply(c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),  function(i)  outtests(Anova(PpoMod1, idata = idata, idesign = ~ block, test.statistic=i)))
tab <- do.call(rbind, tab)

Fval <- tab[,3]
Pval <- tab[,6]
PPOtrue <- data.frame(Fval[2:4], Pval[2:4])
rownames(PPOtrue) <- c("Range", "block", "Range:block")
names(PPOtrue)[names(PPOtrue) == 'Fval.2.4.'] <- 'test.statistic'
names(PPOtrue)[names(PPOtrue) == 'Pval.2.4.'] <- 'p.value'

for (i in 1:3) {
  PPOtrue$BootF[i] <-mean(abs(PpoMat[,i] - mean(PpoMat[,i])) > PPOtrue$test.statistic[i])
}      

# PROTEIN
outtests <- car:::print.Anova.mlm
body(outtests)[[16]] <- quote(invisible(tests))
body(outtests)[[15]] <- NULL
ProMod1 <- lm(cbind(ProteinHr0, ProteinHr4, ProteinHr24) ~ Range, data = ProWide)
tab <- lapply(c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),  function(i)  outtests(Anova(ProMod1, idata = idata, idesign = ~ block, test.statistic=i)))
tab <- do.call(rbind, tab)

Fval <- tab[,3]
Pval <- tab[,6]
PROtrue <- data.frame(Fval[2:4], Pval[2:4])
rownames(PROtrue) <- c("Range", "block", "Range:block")
names(PROtrue)[names(PROtrue) == 'Fval.2.4.'] <- 'test.statistic'
names(PROtrue)[names(PROtrue) == 'Pval.2.4.'] <- 'p.value'

for (i in 1:3) {
  PROtrue$BootF[i] <-mean(abs(ProMat[,i] - mean(ProMat[,i])) > PROtrue$test.statistic[i])
}      

PODtrue
PPOtrue
PROtrue