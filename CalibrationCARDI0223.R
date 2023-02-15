### 1. Preparation
# a. Package loading 
#specify all the packages used in the chapter and install them if they are not already
library(prospectr)
library(soilspec)
library(signal)
library(plyr)
library(pracma)
library(wavethresh)
library(viridis)
library(clhs)
library(pls)
library(reshape2)
library(ggplot2)

# b. Data loading
setwd("C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/4a.Full Averaged databases")
#On charge les donnees
# #
soilspec <- read.table(file = "CARDIFull.csv",
                         head = T, sep = ",", dec = ".", stringsAsFactors = T)


### 2. Organising the database
#Extracting all spectral data
spec <- soilspec[grep("spc", names(soilspec), value = TRUE)]
#remove the spectra from the current dataframe
soilspec <- soilspec[ , -which(names(soilspec) %in% grep("spc", names(soilspec), value = TRUE))]
#add the spectra to the dataframe
soilspec$spc <- spec
#Remove "spc" in front of absorbance value
names(soilspec$spc) = as.numeric(gsub(pattern = "spc.", replacement = "", x = names(soilspec$spc)))


### 3. Pre-processing

#a. Tranform chemical data through log10
# #
soilspec$X.TOC <- log10(soilspec$X.TOC)
#Transform spectral data

#b. Noise Removal for spectra
# #
#b1. Moving Window Average
#specify the window size
windowMa <- 11
# apply the moving average
soilspec$spcMA <- movav(X = soilspec$spcSnvC, w = windowMa)

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
                               m = 2)

#c. Scatter correction for spectra
# #
#c1. Standard Normal Variate
#function for applying standard normal variate transformation
snvBLC <- function(spectra) {
  spectra <- as.matrix(spectra)
  snvMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
  # apply the standardization to each row
  for (i in 1:nrow(spectra)) {
    snvMat[i, ] <- (spectra[i, ] - mean(spectra[i, ]))/sd(spectra[i,])
  }
  snvMat<- as.data.frame(snvMat)
  colnames(snvMat) <- colnames(spectra)
  return(snvMat)
}
# apply the standard normal variate transformation
soilspec$spcSnvC <- snvBLC(soilspec$spc)

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
soilspec$spcMscC <- mscBLC(soilspec$spc)

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

#d. Derivative transformation for spectra
# -> See Savitzky Glay + add value for m

#e.Centring and standardizing
# centre the spectra wavelengths
soilspec$spcNorm <- scale(soilspec$spc, center = TRUE, scale = FALSE)
# standardize the spectra wavelengths
soilspec$spcSdt <- scale(soilspec$spc, center = TRUE, scale = TRUE)

#f. Spectral or Dimension or Reduction

#f1. Resampling
# define the current resolution
oldWavs <- as.numeric(colnames(soilspec$spc))
# define the new resolution, resample every 8 wavelengths
newWavs <- seq(from = min(oldWavs), to = max(oldWavs), by = 5)
# apply the resampling
soilspec$spcResam <- prospectr::resample(X = soilspec$spc,
                                 wav = oldWavs,
                                 new.wav = newWavs,
                                 interpol = "linear")
# dimension of spectra library
dim(soilspec$spc); dim(soilspec$spcResam)

#f2. Wavelets
# function for the wavelet transform of the spectra
# To be developped

#g. Other specific transformation
#g1. Splice correction
#To be developped
#g2. Continuum removal
#To be developped


#Converting data into a matrix
soilspec$spcSg <- as.matrix(soilspec$spcSg)
#Removing of the different outliers
soilspec2 <- subset(soilspec, UnID != "100_020-030" & UnID != "100_030-050"&
                              UnID != "140_000-010" & UnID != "28_000-010" &
                              UnID != "28_010-020"  & UnID != "28_020-030" &
                              UnID != "37_030-050"  & UnID != "37_050-070" &
                              UnID != "37_070-100"  & UnID != "43_000-010" &
                              UnID != "43_010-020"  & UnID != "68_010-020" &
                              UnID != "68_030-050")
                      
                      
#Creation of the different datasets 
TABES <- subset(soilspec2, Cat1 == "ES")
TABPE <- subset(soilspec2, Cat1 == "PE")
TABGPEb <- subset(soilspec2, Cat2 == "GPE"|Cat2 =="GPEb")
TABMPEb <- subset(soilspec2, Cat2 == "MPE"|Cat3 == "MPEb")
TABC1 <- subset(soilspec2, Cat3 == "1")
TABC2 <- subset(soilspec2, Cat3 == "2")

#x. Plotting the spectra
matplot(x = colnames(TABC2$spc), y = t(TABC2$spc),
        xlab = "Wavenumber /cm-1",
        ylab = "Absorbance",
        type = "l",
        xlim = rev(range(as.numeric(colnames(TABC2$spc)))),
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3),
        main = "SG - TABC2")


### 4. Detection of outliers
#To be developped


### 5. Data splitting
# set the seed for reproducibility
# #
set.seed(19101991)
#Set the sample size
SampleSize <- nrow(soilspec)*0.3

#a. Random sampling splitting method
# id of the rows to be used for calibration
ValId <- sample(1:nrow(soilspec), size = round(0.3*nrow(soilspec)))
# separate the dataset into calibration and validation
SSVal <- soilspec[ValId,]
SSCal <- soilspec[-ValId,]

#b. Kennard-stone splitting method
VALid <- kenStone(spec,
                 k = SampleSize,
                 metric = "euclid",
                 .center = TRUE,
                 .scale = FALSE)
ssVal <- soilspec[CALid$test,]
ssCal <- soilspec[-CALid$test,]

#c. K-means splitting method
VALid <- clhs(x = as.data.frame(spec),
              size = SampleSize,
              iter = 1000,
              simple = FALSE)
ssVal <- soilspec[VALid$index_samples,]
ssCal <- soilspec[-VALid$index_samples,]

#d. Duplex splitting method
VALid <- duplex(soilspec$spc,
                k = SampleSize,
                metric = "euclid",
                .center = TRUE,
                .scale = FALSE)
ssVal <- soilspec[VALid$test,]
ssCal <- soilspec[-VALid$test,]

#e. Conditioned Latin Hypercube Sampling
# To be developped


## Check of the chemical data of each sample
# plot the value of the Total Carbon content for both calibration and validation
par(mfrow=c(1,2))
# calibration
hist(ssCal$X.TOC,
     main = "",
     xlab = "Total carbon")
# validation
hist(ssVal$X.TOC,
     main = "",
     xlab = "Total carbon")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

### 6. Calibration through PLSR
## a. PLSR modeling
# maximum number of components in the PLS model
maxc <- 20
# generate a PLS model based on calibration data
soilCPlsModel <- plsr(X.TOC ~ spcSg,
                      data = ssCal,
                      method = "kernelpls",
                      ncomp = maxc,
                      validation = "LOO")
# this is the plsr function, using cross validation to evaluate the RMSEP


##c. Checking the nb of LV comparing results from CAL and VAL then
#Results-Completing Table all indicators with CAL & VAL results
# number of components to use
nc <- 8
# predict on the calibration dataset
soilCPlsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssCal$spcSg)
# predict on the validation dataset
soilVplsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssVal$spcSg)
Cal <- soilspec::eval(ssCal$X.TOC, soilCPlsPred, obj = "quant")
Val <- soilspec::eval(ssVal$X.TOC, soilVplsPred, obj = "quant")
for (i in 1:20){
  # number of components to use
  nc <- i
  # predict on the calibration dataset
  soilCPlsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssCal$spcSg)
  # accuracy measures for calibration
  Cal[i, ] <- soilspec::eval(ssCal$X.TOC, soilCPlsPred, obj = "quant")
  # predict on the validation dataset
  soilVplsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssVal$spcSg)
  # accuracy measures for validation
  Val[i, ] <- soilspec::eval(ssVal$X.TOC, soilVplsPred, obj = "quant")
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
nc <- 8
# predict on the calibration dataset
soilCPlsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssCal$spcSg)
# predict on the validation dataset
soilVplsPred <- predict(soilCPlsModel, ncomp = nc, newdata = ssVal$spcSg)
par(mfrow = c(1, 2))
# plot calibration
plot(ssCal$X.TOC, soilCPlsPred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Calibration data",
     pch = 16)
abline(0, 1)
# plot validation
plot(ssVal$X.TOC, soilVplsPred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Validation data",
     pch = 16)
abline(0, 1)
# Quality model indicators CAL
soilspec::eval(ssCal$X.TOC, soilCPlsPred, obj = "quant")
# Quality model indicators VAL
soilspec::eval(ssVal$X.TOC, soilVplsPred, obj = "quant")

