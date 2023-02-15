## 1. Intro
# On charge les packages
library(plyr)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)
library(tidyverse)
library(stringr)

# On d?fini le dossier de travail
#CARDI Compare 2022 (CAR - C)
setwd("C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/1.Raw spectra/CARDI Compare 2022 (CAR - C)")

## 2. Creation de la base
# Créer une database ? partir de plusieurs fichiers CSV
dataset <- list.files(full.names=TRUE, pattern=".csv") %>% lapply(read.csv, header=TRUE, sep=",") %>% bind_cols()
# Renommer la colonne "Wavelength qu'on va garder"
names(dataset)[names(dataset) == 'Wavenumber...1'] <- 'Wavelength'
# Virer les autres colonnes Wavelength
dataset <- dataset[ , -which(names(dataset) %in% grep("Wavenumber", names(dataset), value = T))]
# Inversion des lignes et colonnes 
dataset <- as.data.frame(t(dataset))
# Placement des longeurs d'onde en nom de colonnes + supression
Wavenumber <- as.list(dataset[1,])
colnames(dataset) <- Wavenumber
dataset <- dataset[-1, ]
# Renommer les colonnes "intensity" par le nom de fichier
newNames  <- list.files(path=".", full.names=TRUE, pattern=".csv")
rownames(dataset) <- newNames
rownames(dataset) = gsub(pattern = "./", replacement = "", x = rownames(dataset))
rownames(dataset) = gsub(pattern = "_0000.csv", replacement = "", x = rownames(dataset))
rownames(dataset) = gsub(pattern = "_0001.csv", replacement = "", x = rownames(dataset))

## 3. Organisation du tableau
# On s?pare l'identifiant en plusieurs colonnes contenant chacune une information
ID <- as.data.frame(rownames(dataset))
colnames(ID) <- "ID"
ID <- ID %>% separate(ID, c("A", "B", "C","D"), sep = "_")
ID$Project <- substr(ID$A, 1, 4)
ID$Sample <- substr(ID$B, 1, 4)
ID$Replicate <- substr(ID$D,4, 4)
ID$SoilDepthDumm <- substr(ID$D,1, 3)
ID$SoilDepth <- paste(ID$C, ID$SoilDepthDumm, sep ="-")
# Supprimer la partie inutile 
ID <- ID[,c(5:7,9)]
# Ordonner les colonnes
ID <- ID %>% relocate(ID, Project)
# Finalement on r?uni les tables
dataset <- cbind(ID,dataset)

# On ordonne le tableau en fonction de l'identifiant
dataset$Sample <- as.numeric(dataset$Sample)
dataset <- dataset[order(dataset[,2],decreasing=F), ]

# Ecart-type intra r?plicat par individus
UnID <- paste(dataset$Sample,dataset$SoilDepth, sep = "_")
dataset <- cbind(UnID,dataset)
SD <- as.data.frame(apply(dataset[,7:length(dataset)],2,tapply, dataset$UnID, sd))
# Moyenne des ?carts types par individus
ID <- rownames(SD)
SDind <- rowMeans(SD)
SDind <- as.data.frame(cbind(ID,SDind))
SDind$SDind <- as.numeric(SDind$SDind)
hist(SDind$SDind, main = "Average standard deviation for each individual", xlab = "Intensity")
# Visualisation des moyennes d'?cart-type par type de sol
ggplot(SDind, aes(ID, SDind, label = ID))+ geom_point()+geom_text(hjust=0, vjust=0)

# Moyenne des ?carts-types par longueur d'onde
SD <- SD[-75,]
WL <- colnames(SD)
WL <- str_replace_all(WL, "X", "")
SDwl <- colMeans(SD)
SDwl <- as.data.frame(t(rbind(WL, SDwl)))
SDwl$SDwl <- as.numeric(SDwl$SDwl)
SDwl$WL <- as.numeric(SDwl$WL)
hist(SDwl$SDwl, main = "Average standard deviation for each wavelength", xlab = "Intensity")
# Visualisation des moyennes d'?cart-type par wavelength
ggplot(SDwl, aes(x=WL, y=SDwl)) +
  geom_line()+ ggtitle("Average standard deviation depending on intensity")+ ylab("Intensity") + xlab("Wavenumber")
ggplot(SDwl, aes(x=WL, y=SDwl,label = WL)) +
  geom_point()+geom_text()+ ggtitle("Average standard deviation depending on intensity")+ ylab("Intensity") + xlab("Wavenumber")
# Elimination des longueurs d'ondes non souhait?es
dataset <- dataset[,c(1:5,27:905)]

## 2. Moyenne des replicats par spectre
# On pr?pare le tableau r?ceptacle
ID <- dataset[,c(1:3,5)]
ID <- ID %>% distinct(UnID, .keep_all = TRUE)

# Calcul des moyennes pour chaque cat?gorie + Ajout d'une seconde colonne indicatrice pour repr?sentation graphique+ Ajout Num
Mean <- as.data.frame(apply(dataset[,6:884],2,tapply, dataset$UnID, mean))
Mean$UnID <- rownames(Mean)
#On r?uni les deux tableaux
DataMean <- merge(ID,Mean, by = "UnID")
#On ?limine les outliers restants (Ici individus: 132_030-050)
DataMean <- DataMean[-41,]

### We add the chemical database
# On d?fini le dossier de travail
setwd("C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/3a.ChemicDatabase")
# On charge les donn?es
dataset <- read.table(file = "CARDI150223.csv", head = T, sep = ";", dec = ".", stringsAsFactors = T)
## 2. Statistiques
# Moyenne & SD
dataset %>%
  group_by(Soil.Depth) %>%
  get_summary_stats(X.TOC, type = "mean_sd")
# Plots
hist(dataset$X.TOC, main = "%C distribution CARDI", xlab = "%")
ggboxplot(dataset, x = "Soil.Depth", y = "X.TOC", color = "Soil.Depth")

### Combinaison des bases de données
DataFull.1 <- merge(dataset,DataMean, by = "UnID")


### 4. Grapique des spectres pour voir si outliers
# a. Mise en forme
# put the spectra into a single dataframe
DataFull.1 <- DataFull.1[,-c(7:9)]
spec <- DataFull.1[,c(7:885)]
# remove the spectra from the current dataframe
soilspec <- DataFull.1[ , -c(7:885)]
# add the spectra to the dataframe
soilspec$spc <- spec
# b. Les Graphiques des spectres averaged
matplot(x = colnames(soilspec$spc), y = t(soilspec$spc),
        xlab = "Wavenumber (cm-1)",
        ylab = "Intensity (Absorbance)",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))
#Graphique des spectres par profondeur
#Couleur par Profondeur
zfac <- factor(DataFull.1$Cat2)
mescouleurs <- rainbow(length(levels(zfac)))
#Graphique g?n?ral
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
## c. #Graphique Spectre par spectre
for (i in 101:200){
  matplot(x = colnames(soilspec$spc[i,]), y = t(soilspec$spc[i,]),
          main = soilspec$UnID[i],
          xlab = "Wavenumber (cm-1)",
          ylab = "Intensity (Absorbance)",
          type = "l",
          lty = 1,
          col = "black",
          ylim = c(0,2))
}



###PCA on soil samples
# Selection des donn?es actives :
TableT.active <- soilspec$spc
# ACP :
res.pca <- PCA(TableT.active, scale.unit = TRUE, ncp = 8, graph = TRUE)
# Scale.unit: une valeur logique. Si TRUE, les donn?es sont standardis?es/normalis?es avant l'analyse.
# Ncp: nombre de dimensions conserv?es dans les r?sultats finaux.
# Graph: une valeur logique. Si TRUE un graphique est affich?.
# Obtention des r?sultats de l'ACP :
# Graphique des valeurs propres
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
# Placement des individus en fct de leur couleur de cos2
fviz_pca_ind (res.pca, col.ind = "cos2",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE # ?vite le chevauchement de texte
)
# Pareil avec des tailles de points
fviz_pca_ind (res.pca, pointsize = "cos2",
              pointshape = 21, fill = "#E7B800",
              repel = TRUE # ?vite le chevauchement de texte
)
#Placement des individus en fonction du groupe
fviz_pca_ind(res.pca, axes = c(1,2),
             geom.ind = "point", # Montre les points seulement (mais pas le "text")
             col.ind = soilspec$Cat2, # colorer by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Ellipses de concentration
             legend.title = "Groups"
)

### Classification Ascendante Hiérarchique
res.hcpc <- HCPC(res.pca)
Clust <- res.hcpc$data.clust
#Intégration des numéros de cluster
soilspec$Cat3 <- Clust[,880]


#Export du tableau final
write.csv(soilspec, file="C:/Users/titou/OneDrive/Bureau/Cambodge 2021/I.VI Carbone du sol/3.Calibration/2. Data/4a.Full Averaged databases/CARDIFull.csv")

