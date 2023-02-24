rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/5-Stats")

# Importation des jeux de données
#################################

patients <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/0-Metadata/patients.txt", 
                        header = TRUE, stringsAsFactors = FALSE )
files <- list.files("~/ZINARA/Gnomic/1-StageM2/1-Appel_de_variant/4-Fichiers_variants_splitEs/")
list_p <- lapply( 1:nrow(patients), function(i) { unlist( strsplit( patients$IDs[i], "," ) ) } )
ID <- unlist( lapply( 1:nrow(patients), function(i) { unlist( strsplit( patients$IDs[i], "," ) ) } ) )
fichiers_par_echantillon <- lapply( 1:length(ID), function(i) { files[ grep( ID[i], files) ] } )
if (!dir.exists("1-text")) { dir.create("1-text") }
if (!dir.exists("2-images")) { dir.create("2-images") }

# Distribution des mutations après appel de variants (avant filtre)
####################################################################

# On met dans un seul dataframe toutes les mutations
dataframe_all <- data.frame( "CHROM"=NA,"POS"=NA,"ID"=NA,"REF"=NA,"ALT"=NA,"QUAL"=NA,"FILTER"=NA,"INFO"=NA,"FORMAT"=NA,
                                       "Details"=NA, "AF"=NA, "Sample"=NA, "Patient"=NA, stringsAsFactors = FALSE )
dataframe_all <- dataframe_all[0,]

library(vcfR)
for (i in 1:27) {
  # On met dans un seul dataframe les mutations filtrées de l'échantillon i
  dataframe_echantillon_i <- data.frame( "CHROM"=NA,"POS"=NA,"ID"=NA,"REF"=NA,"ALT"=NA,"QUAL"=NA,"FILTER"=NA,"INFO"=NA,"FORMAT"=NA,
                                         "Details"=NA, "AF"=NA, "Sample"=NA, "Patient"=NA, stringsAsFactors = FALSE )
  dataframe_echantillon_i <- dataframe_echantillon_i[0,]
  for (fichier in fichiers_par_echantillon[[i]]) {
    # Importation du jeu de donnée
    print( sprintf( "---------------------------Importation du jeu de données %s --------------------------------------", fichier ) )
    data_ <- read.vcfR( file = sprintf( "../1-Appel_de_variant/4-Fichiers_variants_splitEs/%s", fichier ), verbose = FALSE )
    # Vérification si l'échantillon présente des mutations
    print( "     Vérification si l'échantillon présente des mutations" )
    if ( length( data_@fix ) != 0 ) {
      # On remplit la dataframe pour l'échantillon i
      dataframe_<- data.frame( data_@fix, data_@gt, 
                                    extract.gt(x = data_, element = "AF", as.numeric = TRUE), stringsAsFactors = FALSE )
      colnames(dataframe_) <- c( "CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "Details", "AF" )
      # On crée une colonne Sample et une colonne Patient
      dataframe_["Sample"] <- rep( ID[i], nrow(dataframe_) )
      dataframe_["Patient"] <- rep( patients$numero[grep(ID[i], patients$IDs)], nrow(dataframe_) )
      dataframe_echantillon_i <- rbind(dataframe_echantillon_i, dataframe_, stringsAsFactors=F)
      # On remplit la colonne ID avec CHRM+POS+REF+ALT+Sample+Patient
      print( "On remplit la colonne ID avec CHRM+POS+REF+ALT+Sample+Patient" )
      dataframe_echantillon_i$ID <- paste( dataframe_echantillon_i$CHROM, dataframe_echantillon_i$POS, dataframe_echantillon_i$REF,
                                           dataframe_echantillon_i$ALT, dataframe_echantillon_i$Sample, dataframe_echantillon_i$Patient, sep = "_" )
    }
  }
  # On remplit la dataframe de toutes les mutations réunies
  dataframe_all <- rbind(dataframe_all, dataframe_echantillon_i, stringsAsFactors=F)
}

# On crée une colonne DP
dataframe_all["DP"] <- as.numeric( unlist( lapply( 1:nrow(dataframe_all), function(i) { unlist( strsplit( dataframe_all$Details[i], ":" ) )[3] } ) ) )

# On crée une colonne batch
dataframe_all["Batch"] <- ifelse( (substr(dataframe_all$Sample, 1, 3) == "D21"), "batch2", "batch1" )

#####################################################################################################################################################################################
# On crée des un dataframe pour batch1 et un autre pour batch2
data_batch1 <- subset(dataframe_all, batch=="batch1")
data_batch2 <- subset(dataframe_all, batch=="batch2")

# Distribution des DP
print( "On calcule la distribution des DP" )
png(file = "./2-images/1-Histogramme_de_DP.png", width = 900, height = 800)
hist(x = as.numeric(dataframe_all$DP), breaks = 100, main = "", 
     xlab = "Profondeur de lecture (DP)", ylab = "Fréquence", col = "cadetblue" )
dev.off()

# Distribution de DP entre b1 et b2
library(tidyverse)
library(viridis)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()
png(file = "./2-images/1-Boxplot_de_DP_entre_batch_avec_jitter.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=batch, y=as.numeric(DP), fill=batch)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Profondeur de lecture (DP)")
dev.off()

png(file = "./2-images/1-Boxplot_de_DP_entre_batch.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=batch, y=as.numeric(DP), fill=batch)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Profondeur de lecture (DP)")
dev.off()

#data__ <- dataframe_all[order(as.numeric(dataframe_all$Patient)),]
png(file = "./2-images/1-Violin_Chart_DP_entre_batch.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=batch, y=as.numeric(DP), fill=batch)) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=1, option="A") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Profondeur de lecture (DP)")
dev.off()

# Distribution des DP par patient et par échantillon
print( "On calcule la distribution des DP par patient et par échantillon" )
if (!dir.exists("./2-images/1-DP_par_patient_et_par_echantillon_avant_filtre")) { dir.create("./2-images/1-DP_par_patient_et_par_echantillon_avant_filtre") }
for (i in 1:nrow(patients)) {
  for (sample in unlist( strsplit(patients[i,][[2]], ",") )) {
    patient_i_ech_sample_avant_filtre <- subset( dataframe_all, (Patient==i & Sample==sample) )
    png(file = sprintf( "./2-images/1-DP_par_patient_et_par_echantillon_avant_filtre/Histogramme_de_DP_du_patient_%s_echantillon_%s.png", i, sample ), width = 900, height = 800)
    hist(x = as.numeric(patient_i_ech_sample_avant_filtre$DP), breaks = 100, main = sprintf("Patient %s - Echantillon %s", i, sample),
         xlab = "Profondeur de lecture (DP)", ylab = "Fréquence", col = "cadetblue")
    dev.off()
  }
}

# Distribution des DP par patient et par échantillon
print( "On calcule la distribution des DP par patient et par échantillon" )
if (!dir.exists("./2-images/1-DP_par_patient_et_par_echantillon_avant_filtre")) { dir.create("./2-images/1-DP_par_patient_et_par_echantillon_avant_filtre") }
png(file = "./2-images/1-DP_par_patient_et_par_echantillon_avant_filtre/Histogramme_de_DP_par_patient_et_par_echantillon.png", width = 1500, height = 1000)
par(mfrow=c(4, 7))
for (i in 1:nrow(patients)) {
  for (sample in unlist( strsplit(patients[i,][[2]], ",") )) {
    patient_i_ech_sample_avant_filtre <- subset( dataframe_all, (Patient==i & Sample==sample) )
    hist(x = as.numeric(patient_i_ech_sample_avant_filtre$DP), breaks = 100, main = sprintf("Patient %s \n Echantillon %s", i, sample),
         xlab = "Profondeur de lecture (DP)", ylab = "Fréquence", col = "cadetblue")
  }
}
par(mfrow=c(1, 1))
dev.off()

# p <- ggplot(data = dataframe_all, mapping = aes(x = as.character(Patient), y = as.numeric(DP)))
# p + geom_boxplot()
# png(file = "./2-images/1-Boxplot_DP_par_patient.png", width = 900, height = 800)
# p + geom_boxplot(outlier.shape=NA) + coord_flip()
# dev.off()

#data__ <- dataframe_all[order(as.numeric(dataframe_all$Patient)),]
png(file = "./2-images/1-Violin_Chart_DP_par_patient.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=as.character(Patient), y=as.numeric(DP), fill=as.character(Patient))) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=1, option="A") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Patient") + ylab("Profondeur de lecture (DP)")
dev.off()

# Distribution des AF selon DP par échantillon
print( "On calcule la distribution des AF selon DP par échantillon" )
samples_ <- c( "D181208","D202062","D181210","D202061","D181211","D202063","D181212","D202064","D181213","D202065","D181215","D202066","D210285","D210327",
               "D210328","D210288","D210289","D210290","D210297","D210294","D210295","D210296","D210332","D210334","D210335","D210338","D210339" )
if (!dir.exists("./2-images/1-AF_selon_DP_par_echantillon_avant_filtre")) { dir.create("./2-images/1-AF_selon_DP_par_echantillon_avant_filtre") }
for (samp in samples_) {
  echantillon_samp <- subset( dataframe_all, Sample==samp )
  png(file = sprintf( "./2-images/1-AF_selon_DP_par_echantillon_avant_filtre/Nuage_de_points_de_AF_DP_pour_%s.png", samp ), width = 900, height = 800)
  plot(as.numeric(echantillon_samp$AF)~as.numeric(echantillon_samp$DP), main = sprintf("Nuage de points de AF en fonction de DP pour %s", samp),
       xlab = "Profondeur de lecture (DP)", ylab="Fréquence allélique (log(AF))", col = "cadetblue", log="y")
  abline(h=0.3, col="red", )
  dev.off()
}

png(file = "./2-images/1-AF_selon_DP_par_echantillon_avant_filtre/Nuage_de_points_de_AF_DP_tous.png", width = 1500, height = 1000)
par(mfrow=c(4, 7))
for (samp in samples_) {
  echantillon_samp <- subset( dataframe_all, Sample==samp )
  plot(as.numeric(echantillon_samp$AF)~as.numeric(echantillon_samp$DP),
       xlab = "Profondeur de lecture (DP)", ylab="Fréquence allélique (AF)", col = "cadetblue", main = sprintf("%s", samp))
  abline(h=0.3, col="red", )
}
par(mfrow=c(1, 1))
dev.off()

# Distribution des AF selon DP par batch
print( "On calcule la distribution des DP selon AF par batch" )
# batch1 <- c( "D181208","D202062","D181210","D202061","D181211","D202063","D181212","D202064","D181213","D202065","D181215","D202066" )
# batch2 <- c( "D210285","D210327","D210328","D210288","D210289","D210290","D210297","D210294","D210295","D210296","D210332","D210334","D210335","D210338","D210339" )
# data_batch1 <- subset( dataframe_all, (Sample %in% batch1) )
# data_batch2 <- subset( dataframe_all, (Sample %in% batch2) )
data_batch1["AF_OK"] <- ifelse( data_batch1$AF <= 0.3, 1, 0 )
data_batch1["DP_OK"] <- ifelse( data_batch1$DP >= 50, 1, 0 )
data_batch1["AF_DP_OK"] <- ifelse( (data_batch1$AF_OK == 1) & (data_batch1$DP_OK == 1), 1, 0 )
data_batch1$AF_DP_OK <- as.factor(data_batch1$AF_DP_OK)
# my_col <- rainbow(length(levels(data_batch1$AF_DP_OK)))
my_col <- c("grey","red")
if (!dir.exists("./2-images/1-AF_selon_DP_par_batch_avant_filtre")) { dir.create("./2-images/1-AF_selon_DP_par_batch_avant_filtre") }
png(file =  "./2-images/1-AF_selon_DP_par_batch_avant_filtre/Nuage_de_points_de_AF_DP_pour_batch1_log.png", width = 900, height = 800)
plot( as.numeric(data_batch1$DP), as.numeric(data_batch1$AF), main = "Batch 1", log = "xy",
      col = my_col[data_batch1$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
axis(2, 0.3)
axis(1, 50)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

png(file =  "./2-images/1-AF_selon_DP_par_batch_avant_filtre/Nuage_de_points_de_AF_DP_pour_batch1.png", width = 900, height = 800)
plot( as.numeric(data_batch1$DP), as.numeric(data_batch1$AF), main = "Batch 1",
      col = my_col[data_batch1$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

png(file =  "./2-images/1-AF_selon_DP_par_batch_avant_filtre/Nuage_de_points_de_AF_DP_pour_batch1_log_y.png", width = 900, height = 800)
plot( as.numeric(data_batch1$DP), as.numeric(data_batch1$AF), main = "Batch 1", log = "y",
      col = my_col[data_batch1$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

data_batch2["AF_OK"] <- ifelse( data_batch2$AF <= 0.3, 1, 0 )
data_batch2["DP_OK"] <- ifelse( data_batch2$DP >= 50, 1, 0 )
data_batch2["AF_DP_OK"] <- ifelse( (data_batch2$AF_OK == 1) & (data_batch2$DP_OK == 1), 1, 0 )
data_batch2$AF_DP_OK <- as.factor(data_batch2$AF_DP_OK)
# my_col <- rainbow(length(levels(data_batch1$AF_DP_OK)))
my_col2 <- c("grey","red")
png(file =  "./2-images/1-AF_selon_DP_par_batch_avant_filtre/Nuage_de_points_de_AF_DP_pour_batch2_log.png", width = 900, height = 800)
plot( as.numeric(data_batch2$DP), as.numeric(data_batch2$AF), main = "Batch 2", log = "xy",
      col = my_col2[data_batch2$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
axis(2, 0.3)
axis(1, 50)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

png(file =  "./2-images/1-AF_selon_DP_par_batch_avant_filtre/Nuage_de_points_de_AF_DP_pour_batch2.png", width = 900, height = 800)
plot( as.numeric(data_batch2$DP), as.numeric(data_batch2$AF), main = "Batch 2",
      col = my_col2[data_batch2$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

png(file =  "./2-images/1-AF_selon_DP_par_batch_avant_filtre/Nuage_de_points_de_AF_DP_pour_batch2_log_y.png", width = 900, height = 800)
plot( as.numeric(data_batch2$DP), as.numeric(data_batch2$AF), main = "Batch 2", log = "y",
      col = my_col2[data_batch2$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

# Comparaison de DP entre batch1 et batch 2
print( "On fait la comparaison de DP entre batch1 et batch2" )
#p1 <- hist(x = as.numeric(data_batch1$DP), breaks = 100)
#p2 <- hist(x = as.numeric(data_batch2$DP), breaks = 100)
#png(file = "./2-images/1-Comparaison_de_DP_entre_batch1_et_batch2_avant_filtre.png", width = 1200, height = 800)
#plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,50000), main = "Comparaison de DP entre batch1 et batch2")
#plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,50000), add=T)
#legend("right", title="Batch", c("Batch1","Batch2"), fill=c( rgb(0,0,1,1/4), rgb(1,0,0,1/4) ), horiz=FALSE, cex=0.8)
#dev.off()

library(ggplot2)
png(file = "./2-images/1-Comparaison_de_DP_entre_batch1_et_batch2_avant_filtre.png", width = 1200, height = 800)
ggplot(dataframe_all, aes(as.numeric(DP), fill = Batch)) +
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity', binwidth = 30) +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Profondeur de lecture (DP)") + ylab("Densité")
dev.off()

# Distribution des AF par patient et par échantillon
print( "On calcule la distribution des AF par patient et par échantillon" )
if (!dir.exists("./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre")) { dir.create("./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre") }
for (i in 1:nrow(patients)) {
  for (sample in unlist( strsplit(patients[i,][[2]], ",") )) {
    patient_i_ech_sample_avant_filtre <- subset( dataframe_all, (Patient==i & Sample==sample) )
    png(file = sprintf( "./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre/Histogramme_de_AF_du_patient_%s_echantillon_%s.png", i, sample ), width = 900, height = 800)
    hist(x = as.numeric(patient_i_ech_sample_avant_filtre$AF), breaks = 100, main = sprintf("Patient %s - Echantillon %s", i, sample),
         xlab = "Profondeur de lecture (AF)", ylab = "Fréquence", col = "cadetblue")
    dev.off()
  }
}

# Distribution des AF par patient et par échantillon
print( "On calcule la distribution des AF par patient et par échantillon" )
png(file = "./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre/Histogramme_de_AF_par_patient_et_par_echantillon.png", width = 1500, height = 1000)
par(mfrow=c(4, 7))
for (i in 1:nrow(patients)) {
  for (sample in unlist( strsplit(patients[i,][[2]], ",") )) {
    patient_i_ech_sample_avant_filtre <- subset( dataframe_all, (Patient==i & Sample==sample) )
    hist(x = as.numeric(patient_i_ech_sample_avant_filtre$AF), breaks = 100, main = sprintf("Patient %s \n Echantillon %s", i, sample),
         xlab = "Profondeur de lecture (AF)", ylab = "Fréquence", col = "cadetblue")
  }
}
par(mfrow=c(1, 1))
dev.off()

png(file = "./2-images/1-Boxplot_de_AF_entre_batch_avec_jitter.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=batch, y=as.numeric(AF), fill=batch)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Fréquence allélique (AF)")
dev.off()

png(file = "./2-images/1-Boxplot_de_AF_entre_batch.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=batch, y=as.numeric(AF), fill=batch)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Fréquence allélique (AF)")
dev.off()

png(file = "./2-images/1-Violin_Chart_AF_entre_batch.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=batch, y=as.numeric(AF), fill=batch)) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=1, option="A") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Fréquence allélique (AF)")
dev.off()

# Distribution des AF par patient et par échantillon
print( "On calcule la distribution des AF par patient et par échantillon" )
if (!dir.exists("./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre")) { dir.create("./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre") }
for (i in 1:nrow(patients)) {
  for (sample in unlist( strsplit(patients[i,][[2]], ",") )) {
    patient_i_ech_sample_avant_filtre <- subset( dataframe_all, (Patient==i & Sample==sample) )
    png(file = sprintf( "./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre/Histogramme_de_AF_du_patient_%s_echantillon_%s.png", i, sample ), width = 900, height = 800)
    hist(x = as.numeric(patient_i_ech_sample_avant_filtre$AF), breaks = 100, main = sprintf("Patient %s - Echantillon %s", i, sample),
         xlab = "Fréquence allélique (AF)", ylab = "Fréquence", col = "cadetblue")
    dev.off()
  }
}

# Distribution des AF par patient et par échantillon
print( "On calcule la distribution des AF par patient et par échantillon" )
if (!dir.exists("./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre")) { dir.create("./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre") }
png(file = "./2-images/1-AF_par_patient_et_par_echantillon_avant_filtre/Histogramme_de_AF_par_patient_et_par_echantillon.png", width = 1500, height = 1000)
par(mfrow=c(4, 7))
for (i in 1:nrow(patients)) {
  for (sample in unlist( strsplit(patients[i,][[2]], ",") )) {
    patient_i_ech_sample_avant_filtre <- subset( dataframe_all, (Patient==i & Sample==sample) )
    hist(x = as.numeric(patient_i_ech_sample_avant_filtre$AF), breaks = 100, main = sprintf("Patient %s \n Echantillon %s", i, sample),
         xlab = "Fréquence allélique (AF)", ylab = "Fréquence", col = "cadetblue")
  }
}
par(mfrow=c(1, 1))
dev.off()

png(file = "./2-images/1-Violin_Chart_AF_par_patient.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=as.character(Patient), y=as.numeric(AF), fill=as.character(Patient))) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=1, option="A") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Patient") + ylab("Fréquence allélique (AF)")
dev.off()

# Distribution de QUAL
print( "On calcule la distribution de QUAL" )
if (!dir.exists("./2-images/1-QUAL_avant_filtre")) { dir.create("./2-images/1-QUAL_avant_filtre") }
png(file = "./2-images/1-QUAL_avant_filtre/1-Histogramme_de_QUAL_pour_Batch1_et_Batch2.png", width = 1200, height = 1000)
hist( as.numeric(dataframe_all$QUAL), breaks = 100, main = sprintf( "N = %s" , nrow(dataframe_all) ) )
dev.off()
png(file = "./2-images/1-QUAL_avant_filtre/2-Histogramme_de_QUAL_pour_Batch1.png", width = 1200, height = 1000)
hist( as.numeric(data_batch1$QUAL), breaks = 100, main = sprintf( "N = %s ; Pourcentage = %s/100", nrow(data_batch1), round(nrow(data_batch1)/nrow(dataframe_all)*100, 2) ), 
      xlab = "" )
dev.off()
png(file = "./2-images/1-QUAL_avant_filtre/3-Histogramme_de_QUAL_pour_Batch2.png", width = 1200, height = 1000)
hist( as.numeric(data_batch2$QUAL), breaks = 100, main = sprintf( "N = %s ; Pourcentage = %s/100", nrow(data_batch2), round(nrow(data_batch2)/nrow(dataframe_all)*100, 2) ) )
dev.off()

png(file = "./2-images/1-Comparaison_de_QUAL_entre_batch1_et_batch2_avant_filtre.png", width = 1200, height = 800)
ggplot(dataframe_all, aes(as.numeric(QUAL), fill = Batch)) +
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity', binwidth = 30) +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Qualité sur le variant (QUAL)") + ylab("Densité")
dev.off()





















# Distribution des mutations par patient
print( "On calcule la distribution des mutations par patient" )
mutation_per_patient <- table( dataframe_all$Patient )
mutation_per_patient_decreasing = sort( mutation_per_patient, decreasing=TRUE )
mutation_per_patient_decreasing_percentage <- prop.table( mutation_per_patient_decreasing )

png(file = "./2-images/2-Avant_Filtre_par_patient.png", width = 900, height = 800)
bar_plot <- barplot( mutation_per_patient_decreasing_percentage, xlab = "Patient", ylab = "%", ylim = c(0, 0.5), plot = FALSE )
colnames( bar_plot ) <- "x"
barplot( mutation_per_patient_decreasing_percentage, xlab = "Patient", ylab = "%", ylim = c(0, 0.5) )
par(xpd=T)
text(cbind(bar_plot, mutation_per_patient_decreasing_percentage), labels=paste("n=", mutation_per_patient_decreasing, sep = ""), pos=3, offset=0.2, cex = 0.7)
dev.off()

mutation_per_patient_decreasing_percentage_df <- data.frame(mutation_per_patient_decreasing_percentage)
mutation_per_patient_decreasing_df <- data.frame(mutation_per_patient_decreasing)
df1 <- cbind(mutation_per_patient_decreasing_percentage_df, mutation_per_patient_decreasing_df$Freq)
colnames(df1) <- c("Patient", "Pourcentage", "Fréquence")
#df1 <- df1[order(df1$Fréquence),]

df1["Batch"] <- ifelse(as.numeric(levels(df1$Patient)) <= 6, "batch1", "batch2")
png(file = "./2-images/2-Avant_Filtre_par_patient_ggplot.png", width = 900, height = 800)
ggplot(data=df1, aes(x=Patient, y=Pourcentage, fill=Batch)) +
  geom_bar(stat="identity")+
  theme_ipsum() +
  geom_text(aes(label=Fréquence), vjust=-1, color="Black", size=3.5)+
  theme_minimal()
dev.off()

# Distribution des mutations par Batch
print( "On calcule la distribution des mutations par batch" )
mutation_per_batch <- table( dataframe_all$Batch )
mutation_per_batch_decreasing = sort( mutation_per_batch, decreasing = TRUE )
mutation_per_batch_decreasing_percentage <- prop.table( mutation_per_batch_decreasing )

png(file = "./2-images/3-Avant_filtre_par_batch.png", width = 900, height = 800)
bar_plot2 <- barplot( mutation_per_batch_decreasing_percentage, xlab = "Batch", ylab = "%", ylim = c(0, 1), plot = FALSE )
colnames( bar_plot2 ) <- "x"
barplot( mutation_per_batch_decreasing_percentage, xlab = "Batch", ylab = "Pourcentage", ylim = c(0, 1) )
par(xpd=T)
text( cbind( bar_plot2, mutation_per_batch_decreasing_percentage ), labels=paste( "n=", mutation_per_batch_decreasing, sep = "" ), pos=3, offset=0.2, cex=1 )
dev.off()


mutation_per_batch_decreasing_percentage_df <- data.frame(mutation_per_batch_decreasing_percentage)
mutation_per_batch_decreasing_df <- data.frame(mutation_per_batch_decreasing)
df2 <- cbind(mutation_per_batch_decreasing_percentage_df, mutation_per_batch_decreasing_df$Freq)
colnames(df2) <- c("Batch", "Pourcentage", "Fréquence")
png(file = "./2-images/3-Avant_filtre_par_batch_ggplot.pngg", width = 900, height = 800)
ggplot(data=df2, aes(x=Batch, y=Pourcentage, fill=Batch)) +
  geom_bar(stat="identity")+
  theme_ipsum() +
  geom_text(aes(label=Fréquence), vjust=1.6, color="White", size=3.5)+
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
dev.off()


























# Distribution des mutations par échantillon
print( "On calcule la distribution des mutations par échantillon" )
mutation_per_sample <- table(dataframe_all$Sample)
mutation_per_sample_decreasing = sort( mutation_per_sample, decreasing = TRUE )
mutation_per_sample_decreasing_percentage <- prop.table( mutation_per_sample_decreasing )

png(file = "./2-images/4-Avant_filtre_par_echantillon.png", width = 1200, height = 800)
bar_plot3 <- barplot( mutation_per_sample_decreasing_percentage, xlab = "Sample", ylab = "%", plot = FALSE )
colnames( bar_plot3 ) <- "x"
barplot( mutation_per_sample_decreasing_percentage, xlab = "", ylab = "%", las = 2, ylim = c(0,0.25) )
par(xpd=T)
text( cbind( bar_plot3, mutation_per_sample_decreasing_percentage ), labels=paste( "n=", mutation_per_sample_decreasing, sep = "" ), pos=3, offset=0.2, cex=0.8 )
dev.off()

# Distribution des mutations par chromosome
print( "On calcule la distribution des mutations par chromosome" )
mutation_per_chr <- table( dataframe_all$CHROM )
mutation_per_chr_decreasing = sort( mutation_per_chr, decreasing = TRUE )
mutation_per_chr_decreasing_percentage <- prop.table( mutation_per_chr_decreasing )

png( file = "./2-images/5-Avant_filtre_par_chr.png", width = 1200, height = 800 )
bar_plot4 <- barplot( mutation_per_chr_decreasing_percentage, xlab = "Chromosome", ylab = "%", plot = FALSE )
colnames( bar_plot4 ) <- "x"
barplot( mutation_per_chr_decreasing_percentage, xlab = "", ylab = "%", las = 2, ylim = c(0,0.2) )
par(xpd=T)
text( cbind( bar_plot4, mutation_per_chr_decreasing_percentage ), labels=paste( "n=", mutation_per_chr_decreasing, sep = "" ), pos=3, offset=0.2, cex=0.8 )
title(xlab="Chromosome", line=3, cex.lab=1 )
dev.off()



#
write.table(dataframe_all, file="./1-text/Avant_filtre_avec_toutes_les_colonnes.txt", col.names=T, quote=F, row.names=F, sep = "\t")


# Distribution des mutations après filtre
###########################################

# Importation des jeux de données
dataframe_all_patient_prefiltered <- read.table( file = "~/ZINARA/Gnomic/Avec_needlestack_2/Filtre_par_patient_2/Fichiers_filtres/patient_1_filtre.vcf", 
                                                 header = FALSE, stringsAsFactors = FALSE  )
# On crée une colonne sur AF, sur DP, une colonne sur les échantillons, une colonne sur les patients et une colonne sur les batchs
dataframe_all_patient_prefiltered["V11"] <- unlist( lapply( 1:nrow(dataframe_all_patient_prefiltered), function(i) { as.numeric( unlist( strsplit( dataframe_all_patient_prefiltered$V10[i], ":" ) )[6] ) } ) )
dataframe_all_patient_prefiltered["V12"] <- as.numeric( unlist( lapply( 1:nrow(dataframe_all_patient_prefiltered), function(i) { unlist( strsplit( dataframe_all_patient_prefiltered$V10[i], ":" ) )[3] } ) ) )
dataframe_all_patient_prefiltered["V13"] <- unlist( lapply( 1:nrow(dataframe_all_patient_prefiltered), function(i) { unlist( strsplit( dataframe_all_patient_prefiltered$V3[i], "_" ) )[5] } ) )
dataframe_all_patient_prefiltered["V14"] <- rep( "1", nrow(dataframe_all_patient_prefiltered) )
dataframe_all_patient_prefiltered["V15"] <- ifelse( (substr(dataframe_all_patient_prefiltered$V13, 1, 3) == "D21"), "batch2", "batch1" )

for (i in 2:15) {
  
  # On importe le jeu de données suivant à ajouter au précédent
  dataframe_patient_i_prefiltered <- read.table( file = sprintf( "~/ZINARA/Gnomic/Avec_needlestack_2/Filtre_par_patient_2/Fichiers_filtres/patient_%s_filtre.vcf", i ), 
                                                 header = FALSE, stringsAsFactors = FALSE  )
  
  # On crée une colonne sur les patients, une colonne sur les échantillons et une colonne sur les batchs
  dataframe_patient_i_prefiltered["V11"] <- unlist( lapply( 1:nrow(dataframe_patient_i_prefiltered), function(i) { as.numeric( unlist( strsplit( dataframe_patient_i_prefiltered$V10[i], ":" ) )[6] ) } ) )
  dataframe_patient_i_prefiltered["V12"] <- as.numeric( unlist( lapply( 1:nrow(dataframe_patient_i_prefiltered), function(i) { unlist( strsplit( dataframe_patient_i_prefiltered$V10[i], ":" ) )[3] } ) ) )
  dataframe_patient_i_prefiltered["V13"] <- unlist( lapply( 1:nrow(dataframe_patient_i_prefiltered), function(i) { unlist( strsplit( dataframe_patient_i_prefiltered$V3[i], "_" ) )[5] } ) )
  dataframe_patient_i_prefiltered["V14"] <- rep( unlist( strsplit( dataframe_patient_i_prefiltered$V3, "_" ) )[6], nrow(dataframe_patient_i_prefiltered) )
  dataframe_patient_i_prefiltered["V15"] <- ifelse( (substr(dataframe_patient_i_prefiltered$V13, 1, 3) == "D21"), "batch2", "batch1" )
  
  # On combine ce jeu de données au précédent
  dataframe_all_patient_prefiltered <- rbind( dataframe_all_patient_prefiltered, dataframe_patient_i_prefiltered, stringsAsFactors = FALSE )
  
}

# On définit les noms des colonnes
colnames( dataframe_all_patient_prefiltered ) <- c( "CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "Details", "AF", "DP", "Sample", "Patient", "Batch" )

# Distribution des AF par patient et par échantillon
print( "On calcule la distribution des AF par patient et par échantillon" )
if (!dir.exists("./2-images/6-AF_par_patient_et_par_echantillon_apres_filtre")) { dir.create("./2-images/6-AF_par_patient_et_par_echantillon_apres_filtre") }
for (i in 1:nrow(patients)) {
  for (sample in unlist( strsplit(patients[i,][[2]], ",") )) {
    patient_i_ech_sample_apres_filtre <- subset( dataframe_all_patient_prefiltered, (Patient==i & Sample==sample) )
    png(file = sprintf( "./2-images/6-AF_par_patient_et_par_echantillon_apres_filtre/Histogramme_de_AF_du_patient_%s_echantillon_%s.png", i, sample ), width = 900, height = 800)
    hist(x = as.numeric(patient_i_ech_sample_apres_filtre$AF), breaks = 100, main = sprintf("Histogramme de AF du patient %s - échantillon %s", i, sample))
    dev.off()
  }
}

# Distribution des DP par patient et par échantillon
print( "On calcule la distribution des DP par patient et par échantillon" )
if (!dir.exists("./2-images/6-DP_par_patient_et_par_echantillon_apres_filtre")) { dir.create("./2-images/6-DP_par_patient_et_par_echantillon_apres_filtre") }
for (i in 1:nrow(patients)) {
  for (sample in unlist( strsplit(patients[i,][[2]], ",") )) {
    patient_i_ech_sample_apres_filtre <- subset( dataframe_all_patient_prefiltered, (Patient==i & Sample==sample) )
    png(file = sprintf( "./2-images/6-DP_par_patient_et_par_echantillon_apres_filtre/Histogramme_de_DP_du_patient_%s_echantillon_%s.png", i, sample ), width = 900, height = 800)
    hist(x = as.numeric(patient_i_ech_sample_apres_filtre$DP), breaks = 100, main = sprintf("Histogramme de DP du patient %s - échantillon %s", i, sample))
    dev.off()
  }
}

# Distribution des DP
print( "On calcule la distribution des DP" )
png(file = "./2-images/6-Histogramme_de_DP.png", width = 900, height = 800)
hist(x = as.numeric(dataframe_all_patient_prefiltered$DP), breaks = 100, main = "Histogramme de DP" )
dev.off()

# Distribution des AF selon DP par échantillon
print( "On calcule la distribution des AF selon DP par échantillon" )
if (!dir.exists("./2-images/6-AF_selon_DP_par_echantillon_apres_filtre")) { dir.create("./2-images/6-AF_selon_DP_par_echantillon_apres_filtre") }
for (samp in samples_) {
  echantillon_samp <- subset( dataframe_all_patient_prefiltered, Sample==samp )
  png(file = sprintf( "./2-images/6-AF_selon_DP_par_echantillon_apres_filtre/Nuage_de_points_de_AF_DP_pour_%s.png", samp ), width = 900, height = 800)
  plot(as.numeric(echantillon_samp$DP), as.numeric(echantillon_samp$AF), main = sprintf("Nuage de points de AF en fonction de DP pour %s", samp), log = "y")
  dev.off()
}

# Distribution des AF selon DP par batch
print( "On calcule la distribution des AF selon DP par batch" )
data_batch1 <- subset( dataframe_all_patient_prefiltered, Batch == "batch1" )
data_batch2 <- subset( dataframe_all_patient_prefiltered, Batch == "batch2" )
if (!dir.exists("./2-images/6-AF_selon_DP_par_batch_apres_filtre")) { dir.create("./2-images/6-AF_selon_DP_par_batch_apres_filtre") }
png(file =  "./2-images/6-AF_selon_DP_par_batch_apres_filtre/Nuage_de_points_de_AF_DP_pour_batch1.png", width = 900, height = 800)
plot( as.numeric(data_batch1$DP), as.numeric(data_batch1$AF), main = "Nuage de points de AF en fonction de DP pour batch1", log = "y" )
dev.off()
png(file =  "./2-images/6-AF_selon_DP_par_batch_apres_filtre/Nuage_de_points_de_AF_DP_pour_batch2.png", width = 900, height = 800)
plot( as.numeric(data_batch2$DP), as.numeric(data_batch2$AF), main = "Nuage de points de AF en fonction de DP pour batch2", log = "y" )
dev.off()

# Comparaison de DP entre batch1 et batch 2
print( "On fait la comparaison de DP entre batch1 et batch 2" )
#p1 <- hist(x = as.numeric(data_batch1$DP), breaks = 100)
#p2 <- hist(x = as.numeric(data_batch2$DP), breaks = 100)
#png(file = "./2-images/6-Comparaison_de_DP_entre_batch1_et_batch2_apres_filtre.png", width = 1200, height = 800)
#plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,50000), main = "Comparaison de DP entre batch1 et batch2")
#plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,50000), add=T)
#legend("right", title="Batch", c("Batch1","Batch2"), fill=c( rgb(0,0,1,1/4), rgb(1,0,0,1/4) ), horiz=FALSE, cex=0.8)
#dev.off()

dataframe_all_patient_prefiltered["Batch"] <- ifelse( (substr(dataframe_all_patient_prefiltered$Sample, 1, 3) == "D21"), "batch2", "batch1" )
png(file = "./2-images/6-Comparaison_de_DP_entre_batch1_et_batch2_apres_filtre.png", width = 1200, height = 800)
ggplot(dataframe_all_patient_prefiltered, aes(as.numeric(DP), fill = Batch)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity', binwidth = 30)
dev.off()

# Distribution des mutations par patient
print( "On calcule la distribution des mutations par patient" )
mutation_prefiltered_per_patient <- table( dataframe_all_patient_prefiltered$Patient )
mutation_prefiltered_per_patient_decreasing = sort( mutation_prefiltered_per_patient, decreasing=TRUE )
mutation_prefiltered_per_patient_decreasing_percentage <- prop.table( mutation_prefiltered_per_patient_decreasing )

png(file = "./2-images/7-Apres_Filtre_par_patient.png", width = 900, height = 800)
bar_plot5 <- barplot( mutation_prefiltered_per_patient_decreasing_percentage, xlab = "Patient", ylab = "%", ylim = c(0, 0.5), plot = FALSE )
colnames( bar_plot5 ) <- "x"
barplot( mutation_prefiltered_per_patient_decreasing_percentage, xlab = "Patient", ylab = "%", ylim = c(0, 0.5) )
par(xpd=T)
text(cbind(bar_plot5, mutation_prefiltered_per_patient_decreasing_percentage), labels=paste("n=", mutation_prefiltered_per_patient_decreasing, sep = ""), pos=3, offset=0.2, cex = 0.7)
dev.off()

# Distribution des mutations par batch
print( "On calcule la distribution des mutations par batch" )
mutation_prefiltered_per_batch <- table( dataframe_all_patient_prefiltered$Batch )
mutation_prefiltered_per_batch_decreasing = sort( mutation_prefiltered_per_batch, decreasing = TRUE )
mutation_prefiltered_per_batch_decreasing_percentage <- prop.table( mutation_prefiltered_per_batch_decreasing )

png(file = "./2-images/8-Apres_filtre_par_batch.png", width = 900, height = 800)
bar_plot6 <- barplot( mutation_prefiltered_per_batch_decreasing_percentage, xlab = "Batch", ylab = "%", ylim = c(0, 1), plot = FALSE )
colnames( bar_plot6 ) <- "x"
barplot( mutation_prefiltered_per_batch_decreasing_percentage, xlab = "Batch", ylab = "%", ylim = c(0, 1) )
par(xpd=T)
text( cbind( bar_plot6, mutation_prefiltered_per_batch_decreasing_percentage ), labels=paste( "n=", mutation_prefiltered_per_batch_decreasing, sep = "" ), pos=3, offset=0.2, cex=1 )
dev.off()

# Distribution des mutations par échantillon
print( "On calcule la distribution des mutations par échantillon" )
mutation_prefiltered_per_sample <- table(dataframe_all_patient_prefiltered$Sample)
mutation_prefiltered_per_sample_decreasing = sort( mutation_prefiltered_per_sample, decreasing = TRUE )
mutation_prefiltered_per_sample_decreasing_percentage <- prop.table( mutation_prefiltered_per_sample_decreasing )

png(file = "./2-images/9-Apres_filtre_par_echantillon.png", width = 1200, height = 800)
bar_plot7 <- barplot( mutation_prefiltered_per_sample_decreasing_percentage, xlab = "Sample", ylab = "%", plot = FALSE )
colnames( bar_plot7 ) <- "x"
barplot( mutation_prefiltered_per_sample_decreasing_percentage, xlab = "", ylab = "%", las = 2, ylim = c(0,0.25) )
par(xpd=T)
text( cbind( bar_plot7, mutation_prefiltered_per_sample_decreasing_percentage ), labels=paste( "n=", mutation_prefiltered_per_sample_decreasing, sep = "" ), pos=3, offset=0.2, cex=0.8 )
#title(xlab="Sample", line=4, cex.lab=1 )
dev.off()

# Distribution des mutations par chromosome
print( "On calcule la distribution des mutations par chromosome" )
mutation_prefiltered_per_chr <- table( dataframe_all_patient_prefiltered$CHROM )
mutation_prefiltered_per_chr_decreasing = sort( mutation_prefiltered_per_chr, decreasing = TRUE )
mutation_prefiltered_per_chr_decreasing_percentage <- prop.table( mutation_prefiltered_per_chr_decreasing )

png( file = "./2-images/10-Apres_filtre_par_chr.png", width = 1200, height = 800 )
bar_plot8 <- barplot( mutation_prefiltered_per_chr_decreasing_percentage, xlab = "Chromosome", ylab = "%", plot = FALSE )
colnames( bar_plot8 ) <- "x"
barplot( mutation_prefiltered_per_chr_decreasing_percentage, xlab = "", ylab = "%", las = 2, ylim = c(0,0.2) )
par(xpd=T)
text( cbind( bar_plot8, mutation_prefiltered_per_chr_decreasing_percentage ), labels=paste( "n=", mutation_prefiltered_per_chr_decreasing, sep = "" ), pos=3, offset=0.2, cex=0.8 )
title(xlab="Chromosome", line=3, cex.lab=1 )
dev.off()

# Distribution de QUAL
print( "On calcule la distribution de QUAL" )
if (!dir.exists("./2-images/6-QUAL_apres_filtre")) { dir.create("./2-images/6-QUAL_apres_filtre") }
png(file = "./2-images/6-QUAL_apres_filtre/1-Histogramme_de_QUAL_pour_Batch1_et_Batch2.png", width = 1200, height = 1000)
hist( as.numeric(dataframe_all_patient_prefiltered$QUAL), breaks = 100, main = sprintf( "N = %s" , nrow(dataframe_all_patient_prefiltered) ) )
dev.off()
png(file = "./2-images/6-QUAL_apres_filtre/2-Histogramme_de_QUAL_pour_Batch1.png", width = 1200, height = 1000)
hist( as.numeric(data_batch1$QUAL), breaks = 100, main = sprintf( "N = %s ; Pourcentage = %s/100", nrow(data_batch1), round(nrow(data_batch1)/nrow(dataframe_all_patient_prefiltered)*100, 2) ) )
dev.off()
png(file = "./2-images/6-QUAL_apres_filtre/3-Histogramme_de_QUAL_pour_Batch2.png", width = 1200, height = 1000)
hist( as.numeric(data_batch2$QUAL), breaks = 100, main = sprintf( "N = %s ; Pourcentage = %s/100", nrow(data_batch2), round(nrow(data_batch2)/nrow(dataframe_all_patient_prefiltered)*100, 2) ) )
dev.off()

png(file = "./2-images/6-Comparaison_de_QUAL_entre_batch1_et_batch2_apres_filtre.png", width = 1200, height = 800)
ggplot(dataframe_all_patient_prefiltered, aes(as.numeric(QUAL), fill = Batch)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity', binwidth = 30)
dev.off()

#
write.table(dataframe_all_patient_prefiltered, file="./1-text/Apres_filtre_avec_toutes_les_colonnes.txt", col.names=T, quote=F, row.names=F, sep = "\t")


# Distribution des mutations après annotation et filtre sur MAX_AF
#################################################################

# Importation des données
Annotated <- read.table( file = "~/ZINARA/Gnomic/Avec_needlestack_2/Annotation/fichiers_annotes_2_filtres/patient_1_annote.txt", stringsAsFactors = FALSE, sep = "" )
for (i in 2:15) {
  Annotated <- rbind( Annotated, read.table( file = sprintf( "~/ZINARA/Gnomic/Avec_needlestack_2/Annotation/fichiers_annotes_2_filtres/patient_%s_annote.txt", i ), stringsAsFactors = FALSE ), stringsAsFactors = FALSE )
}
colnames( Annotated ) <- c( "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", 
                            "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", 
                            "Codons", "Existing_variation", "Extra" )
write.table(Annotated, file="./1-text/Annotation_fusionnee.txt", col.names=T, quote=F, row.names=F, sep = "\t")

# On crée une colonne échantillon, une colonne patient, une colonne batch et une colonne chromosome
Annotated["Sample"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[5] } ) )
Annotated["Patient"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[6] } ) )
Annotated["Batch"] <- ifelse( (substr(Annotated$Sample, 1, 3) == "D21"), "batch2", "batch1" )
Annotated["Chr"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[1] } ) )

# On crée une colonne Synonyme
#synonymous_variant <- c( "synonymous_variant", "start_retained_variant", "stop_retained_variant" )
synonymous_variant <- c("synonymous_variant", "non_coding_transcript_variant", "start_retained_variant", "synonymous_variant", "non_coding_transcript_exon_variant", "stop_retained_variant")
#nonsynonymous_variant <- c("missense_variant", "redundant_inserted_stop_gained", "start_lost", "stop_gained", "stop_lost")
#nonsynonymous_variant <- c( "nonsynonymous_variant", "missense_variant", "conservative_missense_variant", "non_conservative_missense_variant", "rare_amino_acid_variant", "pyrrolysine_loss", "selenocysteine_loss", "redundant_inserted_stop_gained", "start_lost", "stop_gained", "stop_gained_NMD_escaping", "stop_gained_NMD_triggering", "stop_lost" )
nonsynonymous_variant <- c("nonsynonymous_variant", "conservative_missense_variant", "non_conservative_missense_variant", "stop_gained_NMD_escaping", "stop_gained_NMD_triggering", "rare_amino_acid_variant", "pyrrolysine_loss", "selenocysteine_loss", "redundant_inserted_stop_gained", "frameshift_variant", "inframe_deletion", "splice_acceptor_variant", "stop_lost", "incomplete_terminal_codon_variant", "inframe_insertion", "missense_variant", "protein_altering_variant", "splice_donor_variant", "start_lost", "stop_gained")
Annotated["Synonyme"] <- unlist( lapply( 1:nrow(Annotated), function(i) { ifelse( sum( to_vec(for(j in 1:length( unlist( strsplit( Annotated$Consequence[i], "," ) ) )) unlist( strsplit( Annotated$Consequence[i], "," ) )[j] %in% synonymous_variant ) ) >=1, 1, 0 ) } ) )
Annotated["Nonsynonyme"] <- unlist( lapply( 1:nrow(Annotated), function(i) { ifelse( sum( to_vec(for(j in 1:length( unlist( strsplit( Annotated$Consequence[i], "," ) ) )) unlist( strsplit( Annotated$Consequence[i], "," ) )[j] %in% nonsynonymous_variant ) ) >=1, 1, 0 ) } ) )
write.table(Annotated, file="./1-text/Annotation_fusionnee_avec_colonne_Synonyme_et_Nonsynonyme.txt", col.names=T, quote=F, row.names=F, sep = "\t")

consequences <- unique( Annotated$Consequence )
write.table(consequences, file="./1-text/Consequences_unique.txt", col.names=T, quote=F, row.names=F, sep = "\t")

# Distribution des mutations par patient
#Annotated_mutation_per_patient <- table( Annotated$Patient )
#Annotated_mutation_per_patient_decreasing = sort( Annotated_mutation_per_patient, decreasing=TRUE )
#Annotated_mutation_per_patient_decreasing_percentage <- prop.table( Annotated_mutation_per_patient_decreasing )

#png(file = "./2-images/11-Annotated_per_patient.png", width = 900, height = 800)
#bar_plot9 <- barplot( Annotated_mutation_per_patient_decreasing_percentage, xlab = "Patient", ylab = "%", ylim = c(0, 0.5), plot = FALSE )
#colnames( bar_plot9 ) <- "x"
#barplot( Annotated_mutation_per_patient_decreasing_percentage, xlab = "Patient", ylab = "%", ylim = c(0, 0.5) )
#par(xpd=T)
#text(cbind(bar_plot9, Annotated_mutation_per_patient_decreasing_percentage), labels=paste("n=", Annotated_mutation_per_patient_decreasing, sep = ""), pos=3, offset=0.2, cex = 0.7)
#dev.off()

# Distribution des mutations par batch
#Annotated_mutation_per_batch <- table( Annotated$Batch )
#Annotated_mutation_per_batch_decreasing = sort( Annotated_mutation_per_batch, decreasing = TRUE )
#Annotated_mutation_per_batch_decreasing_percentage <- prop.table( Annotated_mutation_per_batch_decreasing )

#png(file = "./2-images/12-Annotated_per_batch.png", width = 900, height = 800)
#bar_plot10 <- barplot( Annotated_mutation_per_batch_decreasing_percentage, xlab = "Batch", ylab = "%", ylim = c(0, 1), plot = FALSE )
#colnames( bar_plot10 ) <- "x"
#barplot( Annotated_mutation_per_batch_decreasing_percentage, xlab = "Batch", ylab = "%", ylim = c(0, 1) )
#par(xpd=T)
#text( cbind( bar_plot10, Annotated_mutation_per_batch_decreasing_percentage ), labels=paste( "n=", Annotated_mutation_per_batch_decreasing, sep = "" ), pos=3, offset=0.2, cex=1 )
#dev.off()

# Distribution des mutations par échantillon
#Annotated_mutation_per_sample <- table(Annotated$Sample)
#Annotated_mutation_per_sample_decreasing = sort( Annotated_mutation_per_sample, decreasing = TRUE )
#Annotated_mutation_per_sample_decreasing_percentage <- prop.table( Annotated_mutation_per_sample_decreasing )

#png(file = "./2-images/13-Annotated_per_sample.png", width = 1300, height = 800)
#bar_plot11 <- barplot( Annotated_mutation_per_sample_decreasing_percentage, xlab = "Sample", ylab = "%", plot = FALSE )
#colnames( bar_plot11 ) <- "x"
#barplot( Annotated_mutation_per_sample_decreasing_percentage, xlab = "", ylab = "%", las = 2, ylim = c(0,0.25) )
#par(xpd=T)
#text( cbind( bar_plot11, Annotated_mutation_per_sample_decreasing_percentage ), labels=paste( "n=", Annotated_mutation_per_sample_decreasing, sep = "" ), pos=3, offset=0.2, cex=0.8 )
#dev.off()

# Distribution des mutations par chromosome
#Annotated_mutation_per_chr <- table( Annotated$Chr )
#Annotated_mutation_per_chr_decreasing = sort( Annotated_mutation_per_chr, decreasing = TRUE )
#Annotated_mutation_per_chr_decreasing_percentage <- prop.table( Annotated_mutation_per_chr_decreasing )

#png( file = "./2-images/14-Annotated_per_chr.png", width = 1200, height = 800 )
#bar_plot12 <- barplot( Annotated_mutation_per_chr_decreasing_percentage, xlab = "Chromosome", ylab = "%", plot = FALSE )
#colnames( bar_plot12 ) <- "x"
#barplot( Annotated_mutation_per_chr_decreasing_percentage, xlab = "", ylab = "%", las = 2, ylim = c(0,0.12) )
#par(xpd=T)
#text( cbind( bar_plot12, Annotated_mutation_per_chr_decreasing_percentage ), labels=paste( "n=", Annotated_mutation_per_chr_decreasing, sep = "" ), pos=3, offset=0.2, cex=0.8 )
#title(xlab="Chromosome", line=3, cex.lab=1 )
#dev.off()

# Distribution des mutations par gène driver
print( "Importation du jeu de données GENES-DRIVER" )
drivers <- read.table( file = "~/ZINARA/Gnomic/Avec_needlestack_2/Stats_2/IntOGen-DriverGenes_HNSC_2022.tsv", header = TRUE, sep = "\t", dec = ".",
                       stringsAsFactors = FALSE )

# On crée une colonne Gène
print( "On crée une colonne Gène" )
Annotated_reduced <- Annotated
library(stringr)
symbol_position <- str_locate( Annotated$Extra, "SYMBOL=" )
Annotated_reduced["Gene_position_start"] <- symbol_position[,2]
Annotated_reduced <- na.omit(Annotated_reduced)
Annotated_reduced <- transform(Annotated_reduced, Gene_position_start = as.numeric(Gene_position_start)+1)
Annotated_reduced["Extra_reduced"] <- unlist( lapply( 1:nrow(Annotated_reduced), function(i) { substr( Annotated_reduced$Extra[i], as.numeric(Annotated_reduced$Gene_position_start[i]),
                                                                                                       nchar(Annotated_reduced$Extra[i])) } ) )
Gene_position_stop <- str_locate( Annotated_reduced$Extra_reduced, ";" )
Annotated_reduced["Gene_position_stop"] <- Gene_position_stop[,2]-2
Annotated_reduced["Gene"] <- unlist( lapply( 1:nrow(Annotated_reduced), function(i) { substr( Annotated_reduced$Extra[i], as.numeric(Annotated_reduced$Gene_position_start[i]), 
                                                                                              as.numeric(Annotated_reduced$Gene_position_start[i])+as.numeric(Annotated_reduced$Gene_position_stop[i]) ) } ) )

# On calcule le nombre de mutations pour chaque gène driver dans nos jeux de données
print( "On calcule le nombre de mutations pour chaque gène driver dans nos jeux de données" )
driver_list <- data.frame( Symbol = drivers$Symbol, All_mutations = rep( NA, length(drivers$Symbol) ), Synonymous = rep( NA, length(drivers$Symbol) )
                           ,Nonsynonymous = rep( NA, length(drivers$Symbol) ), stringsAsFactors = FALSE ) 

for (i in 1:nrow(driver_list)) {
  #subset_i <- Annotated[ grep( drivers$Symbol[i], Annotated$Extra), ]
  subset_i <- subset( Annotated_reduced, Gene==drivers$Symbol[i] )
  if (nrow(subset_i)!=0) {
    driver_list$All_mutations[i] <- nrow(subset_i)
    driver_list$Nonsynonymous[i] <- sum(as.numeric(subset_i$Nonsynonyme))
    driver_list$Synonymous[i] <- sum(as.numeric(subset_i$Synonyme))
  }
  else
  {
    driver_list$All_mutations[i] <- 0
    driver_list$Nonsynonymous[i] <- 0
    driver_list$Synonymous[i] <- 0
  }
}
driver_list_sorted <- driver_list[order(as.numeric(driver_list$All_mutations), decreasing = TRUE), ]
write.table(driver_list_sorted, file="./1-text/Mutations_synonyme_et_nonsynonyme_par_driver.txt", col.names=T, quote=F, row.names=F, sep = "\t")

# Les 500 gènes les plus mutés parmi les gènes exprimés dans les cellules épithéliales de la muqueuse orale
print( "On cherche les 500 gènes les plus mutés parmi les gènes exprimés dans les cellules épithéliales de la muqueuse orale" )

expressed_genes <- read.table( file = "./1-text/genes_exprimes_dans_les_cellulles_epitheliales_de_la_muqueuse_orale.txt", stringsAsFactors = FALSE )
colnames(expressed_genes) <- "Genes"
expressed_genes["Mutations_number"] <- rep("NA", nrow(expressed_genes))
expressed_genes["Nonsynonyme"] <- rep("NA", nrow(expressed_genes))
expressed_genes["Synonyme"] <- rep("NA", nrow(expressed_genes))
for (i in 1:nrow(expressed_genes)) {
  #subset_i <- Annotated[ grep( expressed_genes$Genes[i], Annotated$Extra), ]
  subset_i <- subset( Annotated_reduced, Gene==expressed_genes$Genes[i] )
  if (nrow(subset_i)!=0) {
    expressed_genes$Mutations_number[i] <- nrow(subset_i)
    expressed_genes$Nonsynonyme[i] <- sum(as.numeric(subset_i$Nonsynonyme))
    expressed_genes$Synonyme[i] <- sum(as.numeric(subset_i$Synonyme))
  }
  else
  {
    expressed_genes$Mutations_number[i] <- 0
    expressed_genes$Nonsynonyme[i] <- 0
    expressed_genes$Synonyme[i] <- 0
  }
}
expressed_genes_sorted <- expressed_genes[order(as.numeric(expressed_genes$Mutations_number), decreasing = TRUE), ]
write.table(expressed_genes_sorted, file="./1-text/genes_exprimes_et_nombre_de_mutations.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table( data.frame(expressed_genes_sorted[1:500,], stringsAsFactors = FALSE), 
             file="./1-text/500_genes_les_plus_mutes_parmi_les_genes_exprimes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

# Les 1000 gènes les plus mutés parmi les gènes non exprimés dans les cellules épithéliales de la muqueuse orale
print( "On cherche les 1000 gènes les plus mutés parmi les gènes non exprimés dans les cellules épithéliales de la muqueuse orale" )

unexpressed_genes <- read.table( file = "./1-text/genes_non_exprimes_dans_les_cellulles_epitheliales_de_la_muqueuse_orale.txt", stringsAsFactors = FALSE )
colnames(unexpressed_genes) <- "Genes"
unexpressed_genes["Mutations_number"] <- rep("NA", nrow(unexpressed_genes))
unexpressed_genes["Nonsynonyme"] <- rep("NA", nrow(unexpressed_genes))
unexpressed_genes["Synonyme"] <- rep("NA", nrow(unexpressed_genes))
for (i in 1:nrow(unexpressed_genes)) {
  #subset_i <- Annotated[ grep( unexpressed_genes$Genes[i], Annotated$Extra), ]
  subset_i <- subset( Annotated_reduced, Gene==unexpressed_genes$Genes[i] )
  if (nrow(subset_i)!=0) {
    unexpressed_genes$Mutations_number[i] <- nrow(subset_i)
    unexpressed_genes$Nonsynonyme[i] <- sum(as.numeric(subset_i$Nonsynonyme))
    unexpressed_genes$Synonyme[i] <- sum(as.numeric(subset_i$Synonyme))
  }
  else
  {
    unexpressed_genes$Mutations_number[i] <- 0
    unexpressed_genes$Nonsynonyme[i] <- 0
    unexpressed_genes$Synonyme[i] <- 0
  }
}
unexpressed_genes_sorted <- unexpressed_genes[order(as.numeric(unexpressed_genes$Mutations_number), decreasing = TRUE), ]
write.table(unexpressed_genes_sorted, file="./1-text/genes_non_exprimes_et_nombre_de_mutations.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table( data.frame(unexpressed_genes_sorted[1:1000,], stringsAsFactors = FALSE), 
             file="./1-text/1000_genes_les_plus_mutes_parmi_les_genes_non_exprimes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

#############################################################################################################################################################
