### I. Preparation

# Package loading
library(plyr)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)
library(tidyverse)
library(prospectr)
library(soilspec)
library(signal)
library(pracma)
library(wavethresh)
library(viridis)
library(clhs)
library(pls)
library(reshape2)
library(ggplot2)
# We define the work directory depending on the data we want to calibrate
  #Global Soil Partnership (GSP - C)
  #Costea Veal Kropeu 2021 (VKC - IC)
  #Rottanak Mondul 2020 (RMP - IC)
  #Stung Chinit 2014 (SCV - IC)
  #Veal Kropeu 2018 (VKV - C)
  #CARDI Compare 2022 (CAR - C)
setwd("C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/1.Raw spectra/Global Soil Partnership (GSP - C)")


### II. RAW database

##A. Database creation
# Create a database from several csv files
dataset <- list.files(full.names=TRUE, pattern=".csv") %>% lapply(read.csv, header=TRUE, sep=",") %>% bind_cols()
# Rename  "Wavelength" column which we will keep
names(dataset)[names(dataset) == 'Wavenumber...1'] <- 'Wavelength'
# Remove other Wavelength columns
dataset <- dataset[ , -which(names(dataset) %in% grep("Wavenumber", names(dataset), value = T))]
# Reverse rows and columns 
dataset <- as.data.frame(t(dataset))
# Add wavenumbers as column names
Wavenumber <- as.list(dataset[1,])
colnames(dataset) <- Wavenumber
dataset <- dataset[-1, ]
# Rename columns "intensity" with csv file names
newNames  <- list.files(path=".", full.names=TRUE, pattern=".csv")
rownames(dataset) <- newNames
rownames(dataset) = gsub(pattern = "./", replacement = "", x = rownames(dataset))
rownames(dataset) = gsub(pattern = "_0000.csv", replacement = "", x = rownames(dataset))
rownames(dataset) = gsub(pattern = "_0001.csv", replacement = "", x = rownames(dataset))

##B. Database organization
# We separate the ID information into several columns 
ID <- as.data.frame(rownames(dataset))
colnames(ID) <- "ID"
ID <- ID %>% separate(ID, c("A", "B"), sep = "-")
ID$Project <- substr(ID$A, 1, 3)
ID$ProjectNumber <- substr(ID$A, 4, 6)
# We separate character chains depending if number or letters
gr <- gregexpr("[0-9//.]+" , ID$B)
ID$ID <- sapply(regmatches(ID$B , gr) , as.numeric)
gr2 <- gregexpr("+[abcde]" , ID$B)
ID$Replicate <- sapply(regmatches(ID$B , gr2) , as.character)
# Removing the useless part 
ID <- ID[,c(3:6)]
# Order the columns
ID <- ID %>% relocate(ID, Project)
# Merge the tables
dataset <- cbind(ID,dataset)
# We order the table depending on the id
dataset$ID <- as.numeric(dataset$ID)
dataset$ProjectNumber <- as.numeric(dataset$ProjectNumber)
dataset2 <- dataset[order(dataset[,1],decreasing=F), ]

##C. Export
# We export the table as a csv file
write.csv(dataset2, file="C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/2a.Raw databases/GSPraw.csv")


###III. INTRA-REPLICATES ANALYZES DATABASE
# Set work directory
setwd("C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/2a.Raw databases")

##A. Data preparation and replicates analyzes
#Load the previous file
  #GSPraw
  #VKCraw
  #RMPraw
  #VKVraw
  #SCVraw
dataset <- read.table(file = "GSPraw.csv", head = T, sep = ",", dec = ".", stringsAsFactors = T)
# Intra-replicate standard deviation per individual
SD <- as.data.frame(apply(dataset[,6:length(dataset)],2,tapply, dataset$ID, sd))
# Average of standard-deviation per individual
ID <- rownames(SD)
SDind <- rowMeans(SD)
SDind <- as.data.frame(cbind(ID,SDind))
SDind$SDind <- as.numeric(SDind$SDind)
hist(SDind$SDind, main = "Average standard deviation for each individual", xlab = "Intensity")
# Visualization of intra replicate standard deviation average value  per individual
ggplot(SDind, aes(ID, SDind, label = ID))+ geom_point()+geom_text(hjust=0, vjust=0)
# Average of standard deviation per wavenumber
WL <- colnames(SD)
WL <- str_replace_all(WL, "X", "")
SDwl <- colMeans(SD)
SDwl <- as.data.frame(t(rbind(WL, SDwl)))
SDwl$SDwl <- as.numeric(SDwl$SDwl)
SDwl$WL <- as.numeric(SDwl$WL)
hist(SDwl$SDwl, main = "Average standard deviation for each wavelength", xlab = "Intensity")
# Visualization of average standard deviation per wavelength
ggplot(SDwl, aes(x=WL, y=SDwl)) +
  geom_line()+ ggtitle("Average standard deviation depending on intensity")+ ylab("Intensity") + xlab("Wavenumber")
ggplot(SDwl, aes(x=WL, y=SDwl,label = WL)) +
  geom_point()+geom_text()+ ggtitle("Average standard deviation depending on intensity")+ ylab("Intensity") + xlab("Wavenumber")

##B. Outlier removing and exporing
# Removing of non-desired wavelengths
dataset2 <- dataset[,c(1:5,27:905)]
# Removing of individuals with high intra-replicate variability
dataset2 <- subset(dataset2, ID != 872 & ID != 970 & ID != 1005 & ID != 1006 & ID != 1026)
# Data export as csv file
write.csv(dataset2, file="C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/2b. Raw databases Clean/GSPclean.csv")


### IV. DATA AVERAGING
# Set work directory
setwd("C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/2b. Raw databases Clean")

##A. Calculation of mean value for each spectra
# Importation
  #GSPclean
  #VKCclean
  #RMPclean
  #VKVclean
  #SCVclean
dataset <- read.table(file = "GSPclean.csv", head = T, sep = ",", dec = ".", stringsAsFactors = T)
rownames(dataset) <- dataset$X
# Table preparation
ID <- dataset[,c(1:3)]
ID <- ID %>% distinct(ID, .keep_all = TRUE)
# Calculation of the mean values
Mean <- as.data.frame(apply(dataset[,7:885],2,tapply, dataset$ID, mean))
#Table merging
DataMean <- cbind(ID,Mean)
#We remove "a" letter from columns
rownames(DataMean) <- str_replace_all(rownames(DataMean), "a", "")
DataMean$X <- str_replace_all(rownames(DataMean), "a", "")

##B. We add the table including chemical analyzes values
# On importe le tableau
setwd("C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/3a.ChemicDatabase")
Chemic <- read.csv2(file = "DB12052022.csv", head = T, sep = ";", dec = ".", stringsAsFactors = T)
#Select only the rigth columns and rows
  #DALRM - GSP FAO
  #COSTEA - Veal Kropeu
  #Sangha village / CASF
  #Stung chinit
  #Veal Kropeu village
Chemic <- subset(Chemic, X1.Source..from.which.village.project.station. == "DALRM - GSP FAO")
#Table preparation
DataMean$Jointure <- rownames(DataMean)
DataMean <- DataMean %>% relocate(Jointure, ID)
Chemic$Jointure <- Chemic$Spectra.database.ID
#We remove the "X" symbols
colnames(Chemic) <- str_replace_all(colnames(Chemic), "X", "A")
#We merge the chemical table with spectra values
DataFull.1 <- left_join(Chemic, DataMean, by="Jointure")
DataFull.1 <- subset(DataFull.1, X726.8309 != "")
#We remove the useless columns
DataFull.1 <- DataFull.1[,c(1:69,72,74:953)]


### V. Plotting the spectra

##A. Table preparation
# put the spectra into a single dataframe
dum <- colnames(DataFull.1[,c(70:950)])
dum[1:3] <- c("Name", "Nb", "Name3")
colnames(DataFull.1[,c(70:950)]) <- dum
spec <- DataFull.1[grep("X", names(DataFull.1), value = TRUE)]
spec <- spec[,-c(1:2)]
# remove the spectra from the current dataframe
soilspec <- DataFull.1[ , -which(names(DataFull.1) %in% grep("X", names(DataMean), value = TRUE))]
# add the spectra to the dataframe
soilspec$spc <- spec
#b. Remove X in front of absorbance value
# take each column name from the spectra DataMean
oldNames <- grep("X", names(soilspec$spc), value = TRUE)
# remove the "X" and make a numeric vector
wavelength <- as.numeric(substring(grep("X", names(soilspec$spc), value = T), 2, 20))
# change the name of the columns of the spectra
colnames(soilspec$spc) <- wavelength
# b. Averaged spectra plot
matplot(x = colnames(soilspec$spc), y = t(soilspec$spc), 
        xlab = "Wavenumber (cm-1)",
        ylab = "Intensity (Absorbance)",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))
#Spectra plot depending on soil depth
zfac <- factor(DataFull.1$A2.Depth..cm.)
mescouleurs <- rainbow(length(levels(zfac)))
#General plot
matplot(x = colnames(soilspec$spc), y = t(soilspec$spc),
        xlab = "Wavenumber (cm-1)",
        ylab = "Intensity (Absorbance)",
        type = "l",
        lty = 1,
        col = mescouleurs[zfac], )
legend("top",
       col = mescouleurs,
       pch = 19,
       legend = levels(zfac))

##B. Data export
write.csv(DataFull.1, file="C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/4a.Full Averaged databases/GSPFull.csv")


### VI. CALIBRATION

##A. Preparation
#Set work directory
setwd("C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/4a.Full Averaged databases")
#Data loading
soilspec <- read.table(file = "GSPFull.csv",
                       head = T, sep = ",", dec = ".", stringsAsFactors = T)
#Extracting all spectral data
spec <- soilspec[grep("X", names(soilspec), value = TRUE)]
spec <- spec[,-c(1:2)]
#remove the spectra from the current dataframe
soilspec <- soilspec[ , -which(names(soilspec) %in% grep("X", names(soilspec),
                                                         value = TRUE))]
#add the spectra to the dataframe
soilspec$spc <- spec
#Remove X in front of absorbance value
oldNames <- grep("X", names(soilspec$spc), value = TRUE)
# remove the "X" and make a numeric vector
wavelength <- as.numeric(substring(grep("X", names(soilspec$spc), value = T), 2, 20))
# change the name of the columns of the spectra
colnames(soilspec$spc) <- wavelength


##B.Pre-processing
#a. Tranform chemical data through log10
# #
soilspec$logA3.SOC.w..wet.combustion. <- log10(soilspec$A3.SOC.w..wet.combustion.)
#Transform spectral data


#b2. Savitzky-Golay Filtering
#function for applying Savitzky-Golay smoothing filter 
#w = windows size, k = the polynomial order; m = derivative of the polynomial
filterSg <- function(spectra, w, k, m) {
  spectra <- as.matrix(spectra)
  ## run filter, the window size in the sgolayfilt function is called n and
  ## the polynomial order is called p
  sg <- aaply(spectra, 1, sgolayfilt, n = w, p = k, m = m)
  ## arrange appropriately if a single sample
  if (nrow(spectra) == 1) {
    sg <- matrix(sg, dim(spectra))
  }
  ## return data frame
  sg <- as.data.frame(sg)
  colnames(sg) <- colnames(spectra)
  return(sg)
}
# filter spectra dataset
soilspec$spcSg <- filterSg(soilspec$spc,
                           w = 11,
                           k = 2,
                           m = 1)

#c. Scatter correction for spectra

#c2. Multiplicative Scatter Correction
# function for applying multiplicative scatter correction
mscBLC <- function(spectra) {
  # first calculate a mean spectrum.
  meanSpec <- as.matrix(colMeans(spectra))
  mscMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
  spectra <- as.matrix(spectra)
  # make a loop over each row
  for (i in 1:nrow(spectra)) {# determine the slope and intercept coefficients
    specLM <- lm(spectra[i, ] ~ meanSpec)
    specCE <- t(as.matrix(specLM$coefficients))# adjust the spectra
    mscMat[i, ] <- t(as.matrix((spectra[i, ] - specCE[1, 1])/specCE[1,2]))
  }
  mscMat <- as.data.frame(mscMat)
  colnames(mscMat) <- colnames(spectra)
  return(mscMat)
}
# apply the multiplicative scatter correction
soilspec$spcMsc <- mscBLC(soilspec$spc)

#c3. Detrend 
# function for detrending a matrix of spectra
detrendSpc <- function(spectra) {
  detrendMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
  spectra <- as.matrix(spectra)
  # make a loop over each row
  for (i in 1:nrow(spectra)) {# detrend each spectra, specify the linear model
    specLM <- pracma::detrend(spectra[i, ], tt = "linear")# take the values and store in the matrix
    detrendMat[i, ] <- as.numeric(specLM[,1])
  }
  detrendMat<- as.data.frame(detrendMat)
  colnames(detrendMat) <- colnames(spectra)
  return(detrendMat)
}
# detrend the spectra
soilspec$spcDT <- detrendSpc(soilspec$spc)


#x. Plotting the spectra
# #
matplot(x = colnames(soilspec$spc), y = t(soilspec$spc),
        xlab = "Wavenumber /cm-1",
        ylab = "Absorbance",
        type = "l",
        xlim = rev(range(as.numeric(colnames(soilspec$spcDT)))),
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
        main = "NPT")


### 5. Data splitting
# set the seed for reproducibility
# #
set.seed(19101991)
#Set the sample size
SampleSize <- nrow(soilspec)*0.3


#d. Duplex splitting method
VALid <- duplex(spec,
                k = SampleSize,
                metric = "euclid",
                .center = TRUE,
                .scale = FALSE)
ssVal <- soilspec[VALid$test,]
ssCal <- soilspec[-VALid$test,]


## Check of the chemical data of each sample
# plot the value of the Total Carbon content for both calibration and validation
par(mfrow=c(1,2))
# calibration
hist(ssCal$logA3.SOC.w..wet.combustion.,
     main = "CAL",
     xlab = "Total carbon")
# validation
hist(ssVal$logA3.SOC.w..wet.combustion.,
     main = "VAL",
     xlab = "Total carbon")


### 6. Calibration through PLSR
## a. PLSR modeling
# maximum number of components in the PLS model
ssCal$spcSg <- as.matrix(ssCal$spcSg)
ssVal$spcSg <- as.matrix(ssVal$spcSg)
maxc <- 20
# generate a PLS model based on calibration data
soilCPlsModel <- plsr(logA3.SOC.w..wet.combustion. ~ spcSg,
                      data = ssCal,
                      method = "kernelpls",
                      ncomp = maxc,
                      validation = "LOO")
# this is the plsr function, using cross validation to evaluate the RMSEP

# Choice of nb of LV
plot(soilCPlsModel, "val",
     main = " ",
     xlab = "Number of components")
ncomp.onesigma <- selectNcomp(soilCPlsModel, method = "onesigma", plot = TRUE)
ncomp.permut <- selectNcomp(soilCPlsModel, method = "randomization", plot = TRUE)

##c. Checking the nb of LV comparing results from CAL and VAL then
#Results-Completing Table all indicators with CAL & VAL results
# number of components to use
nc <- 6
# predict on the calibration dataset
soilCPlsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssCal$spcSg)
# predict on the validation dataset
soilVplsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssVal$spcSg)
Cal <- soilspec::eval(ssCal$logA3.SOC.w..wet.combustion., soilCPlsPred, obj = "quant")
Val <- soilspec::eval(ssVal$logA3.SOC.w..wet.combustion., soilVplsPred, obj = "quant")
for (i in 1:20){
  # number of components to use
  nc <- i
  # predict on the calibration dataset
  soilCPlsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssCal$spcSg)
  # accuracy measures for calibration
  Cal[i, ] <- soilspec::eval(ssCal$logA3.SOC.w..wet.combustion., soilCPlsPred, obj = "quant")
  # predict on the validation dataset
  soilVplsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssVal$spcSg)
  # accuracy measures for validation
  Val[i, ] <- soilspec::eval(ssVal$logA3.SOC.w..wet.combustion., soilVplsPred, obj = "quant")
}
#Merge the two tables
R2RMSE <- cbind(Cal,Val)
colnames(R2RMSE) <- c("CAL-ME", "CAL-RMSE", "CAL-r2", "CAL-R2", "CAL-rhoC",
                      "CAL-RPD", "CAL-RPIQ","VAL-ME", "VAL-RMSE", "VAL-r2",
                      "VAL-R2", "VAL-rhoC", "VAL-RPD", "VAL-RPIQ")
#Plots RMSE CAL vs VAL
RMSE <- cbind(R2RMSE$`CAL-RMSE`,R2RMSE$`VAL-RMSE`)
colnames(RMSE) <- c("CAL","VAL")
RMSE.melted <- melt(RMSE, id="x")
qplot(x=Var1, y=value, color=Var2, data=RMSE.melted , geom="line",
      ylab = "RMSE", xlab = "NB LV", main = "RMSE CAL vs VAL")
#Plots R2 CAL vs VAL
R2 <- cbind(R2RMSE$`CAL-R2`,R2RMSE$`VAL-R2`)
colnames(R2) <- c("CAL","VAL")
R2.melted <- melt(R2, id="x")
qplot(x=Var1, y=value, color=Var2, data=R2.melted , geom="line",
      ylab = "RMSE", xlab = "NB LV", main = "RMSE CAL vs VAL")

## d. Plotting obs vs pred for chosen nb of LV, CAL and VAL
# number of components to use
nc <- 4
# predict on the calibration dataset
soilCPlsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssCal$spcSg)
# predict on the validation dataset
soilVplsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssVal$spcSg)
par(mfrow = c(1, 2))
# plot calibration
plot(ssCal$logA3.SOC.w..wet.combustion., soilCPlsPred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Calibration data",
     pch = 16)
abline(0, 1)
# plot validation
plot(ssVal$logA3.SOC.w..wet.combustion., soilVplsPred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Validation data",
     pch = 16)
abline(0, 1)
# Quality model indicators CAL
soilspec::eval(ssCal$logA3.SOC.w..wet.combustion., soilCPlsPred, obj = "quant")
# Quality model indicators VAL
soilspec::eval(ssVal$logA3.SOC.w..wet.combustion., soilVplsPred, obj = "quant")


