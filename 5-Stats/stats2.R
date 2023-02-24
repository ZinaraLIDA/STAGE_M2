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
      dataframe_["Patient"] <- rep( as.character(patients$numero[grep(ID[i], patients$IDs)]), nrow(dataframe_) )
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

# On crée une colonne AO
dataframe_all["AO"] <- unlist(lapply( 1:nrow(dataframe_all), function(i) { as.numeric(unlist(strsplit(dataframe_all$Details[i], ":"))[5]) } ))

#####################################################################################################################################################################################
# On crée des un dataframe pour batch1 et un autre pour batch2
dataframe_all$QUAL <- as.numeric(dataframe_all$QUAL )
data_batch1 <- subset(dataframe_all, Batch=="batch1")
data_batch2 <- subset(dataframe_all, Batch=="batch2")

# Distribution des DP
png(file = "./2-images/1-Histogramme_et_boxplot_de_DP.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 4, 1.1, 2.1))
boxplot(as.numeric(dataframe_all$DP), horizontal=TRUE, ylim=c(0,max(as.numeric(dataframe_all$DP))), xaxt="n",
        col=rgb(0.8,0.8,0,0.5), frame=F, cex = 0.3, pch=19)
par(mar=c(4, 4, 1.1, 2.1))
hist(as.numeric(dataframe_all$DP) , breaks=100, col="cadetblue", border=F, main="" , 
     xlab="Profondeur de lecture (DP)", ylab = "Fréquence", xlim=c(0,max(as.numeric(dataframe_all$DP))))
dev.off()

# Distribution de DP entre b1 et b2
test_ <- wilcox.test(data_batch1$DP, data_batch2$DP, paired = FALSE)
p <- "< 2.2e-16"
png(file = "./2-images/2-Histogramme_et_boxplot_de_DP_entre_batch.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 5.5, 1.1, 2.1))
boxplot(log10(as.numeric(dataframe_all$DP))~dataframe_all$Batch, horizontal=TRUE, xaxt="n",
        col=c("brown1", "#00AFBB"), frame=F, cex = 0.3, pch=19, axes=FALSE, ann = FALSE,
        ylim=c(0,6))
par(mar=c(4, 5.5, 0, 2.1))
options(scipen=100)
hist(log10(as.numeric(data_batch2$DP)), breaks=25, xlim=c(0,6),
     col="#00AFBB", xlab="Profondeur de lecture (DP) en log10", border=F,
     ylab="Fréquence", main="", cex.main = 2,
     cex.lab = 3 , cex.axis=2.5, ylim = c(0, 600000) )
hist(log10(data_batch1$DP), breaks=25, xlim=c(0,6), border=F,
     col="brown1", add=T)
# legend("right", legend=c("Batch2","Batch1"), col=c("#00AFBB", "brown1"),
#        pt.cex=2, pch=15, border = "white", bty = "n", cex = 2)
dev.off()

library(tidyverse)
library(viridis)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()
png(file = "./2-images/3-Violin_Chart_DP_entre_batch.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Batch, y=as.numeric(DP), fill=Batch)) +
  geom_violin() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Profondeur de lecture (DP)")
dev.off()

# Distribution des DP par patient
png(file = "./2-images/4-Violin_Chart_DP_par_patient.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Patient, y=as.numeric(DP), fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Patient") + ylab("Profondeur de lecture (DP)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()

# Distribution des DP par échantillon
png(file = "./2-images/5-Violin_Chart_DP_par_echantillon.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Sample, y=as.numeric(DP), fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(angle = 45)
  ) +
  xlab("Echantillon") + ylab("Profondeur de lecture (DP)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
  #scale_color_viridis(discrete = TRUE, alpha=1, option="A")
  #scale_color_brewer(palette="Dark2")
  # scale_color_manual(values=c("#000000", "#FF0033", "#99CC00", "#FFFF00", "#996600",
  #                             "#FF00CC", "#FF9933", "#0099FF", "#999999", "#990000",
  #                             "#336600", "#CCCC33", "#663300", "#9900CC", "#0033FF"))
dev.off()

####################################################################################################################
# Distribution des AF
png(file = "./2-images/6-Histogramme_et_boxplot_de_AF.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 4, 1.1, 2.1))
boxplot(as.numeric(dataframe_all$AF), horizontal=TRUE, ylim=c(0,max(as.numeric(dataframe_all$AF))), xaxt="n",
        col=rgb(0.8,0.8,0,0.5), frame=F, cex = 0.3, pch=19)
par(mar=c(4, 4, 1.1, 2.1))
hist(as.numeric(dataframe_all$AF) , breaks=100, col="cadetblue", border=F, main="" , 
     xlab="Fréquence allélique (AF)", ylab = "Fréquence", xlim=c(0,max(as.numeric(dataframe_all$AF))))
dev.off()


png(file = "./2-images/6-Histogramme_AF_qlq_exemples.png", width = 1000, height = 800)
par(mfrow=c(2,2))
samples <- c("D181208", "D181210", "D181211", "D181212")
for (samp in samples) {
  pat <- subset( dataframe_all, Sample==samp )
  hist(x = as.numeric(pat$AF), breaks = 100,
       main = sprintf("Histogramme de AF de l'échantillon %s", samp), freq = TRUE, 
       xlab = "Fréquence allélique (AF)", ylab = "Fréquance" )
}
dev.off()
par(mfrow=c(1,1))


# Distribution de AF entre b1 et b2
test_ <- wilcox.test(data_batch1$AF, data_batch2$AF, paired = FALSE)
p <- "< 2.2e-16"
png(file = "./2-images/7-Histogramme_et_boxplot_de_AF_entre_batch.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 5.5, 1.1, 2.1))
boxplot(as.numeric(dataframe_all$AF)~dataframe_all$Batch, horizontal=TRUE, xaxt="n",
        col=c("brown1", "#00AFBB"), frame=F, cex = 0.3, pch=19, axes=FALSE, ann = FALSE)
par(mar=c(4, 5.5, 0, 2.1))
hist(as.numeric(data_batch2$AF), breaks=30, xlim=c(0,max(as.numeric(dataframe_all$AF))),
     col="#00AFBB", xlab="Fréquence allélique (AF)", border=F,
     ylab="Fréquence", cex.lab=3, ylim=c(0, 2000000),
     main="", cex.axis=2.5 )
hist(data_batch1$AF, breaks=30, xlim=c(0,max(as.numeric(dataframe_all$AF))), border=F,
     col="brown1", add=T)
legend("right", legend=c("batch2","batch1"), col=c("#00AFBB", "brown1"),
       pt.cex=5, pch=15, border = "white", bty = "n", cex = 3)
dev.off()

library(tidyverse)
library(viridis)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()
png(file = "./2-images/8-Violin_Chart_AF_entre_batch.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Batch, y=as.numeric(AF), fill=Batch)) +
  geom_violin() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Fréquence allélique (AF)")
dev.off()

# Distribution des AF par patient
png(file = "./2-images/9-Violin_Chart_AF_par_patient.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Patient, y=as.numeric(AF), fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Patient") + ylab("Fréquence allélique (AF)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()

# Distribution des AF par échantillon
png(file = "./2-images/10-Violin_Chart_AF_par_echantillon.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Sample, y=as.numeric(AF), fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(angle = 45)
  ) +
  xlab("Echantillon") + ylab("Fréquence allélique (AF)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
#scale_color_viridis(discrete = TRUE, alpha=1, option="A")
#scale_color_brewer(palette="Dark2")
# scale_color_manual(values=c("#000000", "#FF0033", "#99CC00", "#FFFF00", "#996600",
#                             "#FF00CC", "#FF9933", "#0099FF", "#999999", "#990000",
#                             "#336600", "#CCCC33", "#663300", "#9900CC", "#0033FF"))
dev.off()
####################################################################################################################
# Distribution des AF selon DP par batch
my_col <- c("grey","red")
data_batch1["AF_OK"] <- ifelse( data_batch1$AF <= 0.3, 1, 0 )
data_batch1["DP_OK"] <- ifelse( data_batch1$DP >= 50, 1, 0 )
data_batch1["AF_DP_OK"] <- ifelse( (data_batch1$AF_OK == 1) & (data_batch1$DP_OK == 1), 1, 0 )
data_batch1$AF_DP_OK <- as.factor(data_batch1$AF_DP_OK)
png(file =  "./2-images/11-Nuage_de_points_de_AF_DP_pour_batch1.png", width = 900, height = 800)
plot( as.numeric(data_batch1$DP), as.numeric(data_batch1$AF), main = "Batch 1",
      col = my_col[data_batch1$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (DP)", ylab = "Fréquence allélique (AF)", cex = 0.5)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

data_batch2["AF_OK"] <- ifelse( data_batch2$AF <= 0.3, 1, 0 )
data_batch2["DP_OK"] <- ifelse( data_batch2$DP >= 50, 1, 0 )
data_batch2["AF_DP_OK"] <- ifelse( (data_batch2$AF_OK == 1) & (data_batch2$DP_OK == 1), 1, 0 )
data_batch2$AF_DP_OK <- as.factor(data_batch2$AF_DP_OK)
png(file =  "./2-images/12-Nuage_de_points_de_AF_DP_pour_batch2.png", width = 900, height = 800)
plot( as.numeric(data_batch2$DP), as.numeric(data_batch2$AF), main = "Batch 2",
      col = my_col[data_batch2$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (DP)", ylab = "Fréquence allélique (AF)", cex = 0.5)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

####

png(file =  "./2-images/13-Nuage_de_points_de_AF_DP_pour_batch1_log_y.png", width = 900, height = 800)
plot( as.numeric(data_batch1$DP), as.numeric(data_batch1$AF), main = "Batch 1", log = "y",
      col = my_col[data_batch1$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (DP)", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

png(file =  "./2-images/14-Nuage_de_points_de_AF_DP_pour_batch2_log_y.png", width = 900, height = 800)
plot( as.numeric(data_batch2$DP), as.numeric(data_batch2$AF), main = "Batch 2", log = "y",
      col = my_col[data_batch2$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (DP)", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

####

my_col <- c("grey","red")
png(file =  "./2-images/15-Nuage_de_points_de_AF_DP_pour_batch1_log.png", width = 900, height = 800)
plot( as.numeric(data_batch1$DP), as.numeric(data_batch1$AF), main = "Batch 1", log = "xy",
      col = my_col[data_batch1$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
axis(2, 0.3)
axis(1, 50)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()

png(file =  "./2-images/16-Nuage_de_points_de_AF_DP_pour_batch2_log.png", width = 900, height = 800)
plot( as.numeric(data_batch2$DP), as.numeric(data_batch2$AF), main = "Batch 2", log = "xy",
      col = my_col[data_batch2$AF_DP_OK], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
axis(2, 0.3)
axis(1, 50)
abline(h=0.3, v=50, col="blue", lwd=2)
dev.off()
 ####

dataframe_all["AF_DP_AO_QUAL"] <- ifelse( ((dataframe_all$Batch=="batch1") & (dataframe_all$AF<=0.3) & (dataframe_all$DP>=50) & (dataframe_all$AO>=3) & (dataframe_all$QUAL>=100)) |
                                            ((dataframe_all$Batch=="batch2") & (dataframe_all$AF<=0.3) & (dataframe_all$DP>=50) & (dataframe_all$AO>=3) & (dataframe_all$QUAL>=400)),
                                          1, 0)
dataframe_all$AF_DP_AO_QUAL <- as.factor(dataframe_all$AF_DP_AO_QUAL)

png(file =  "./2-images/16-Nuage_de_points_de_AF_DP_AO_QUAL_log.png", width = 900, height = 800)
par(mar=c(5, 6, 5, 1))
plot( as.numeric(dataframe_all$DP), as.numeric(dataframe_all$AF), main = "AF<0.3 ; DP>50 ; AO>3 \n QUALbatch1>100 ; QUALbatch2>400", log = "xy",
      col = my_col[dataframe_all$AF_DP_AO_QUAL], pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))",
      cex = 0.5, cex.lab=3, cex.axis=2.2, cex.main=2.5)
axis(2, 0.3)
axis(1, 50)
abline(h=0.3, v=50, col="blue", lwd=2, cex=2.5)
dev.off()

####################################################################################################################
# Distribution des QUAL
png(file = "./2-images/17-Histogramme_et_boxplot_de_QUAL.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 4, 1.1, 2.1))
boxplot(dataframe_all$QUAL, horizontal=TRUE, ylim=c(0,max(dataframe_all$QUAL)), xaxt="n",
        col=rgb(0.8,0.8,0,0.5), frame=F, cex = 0.3, pch=19)
par(mar=c(4, 4, 1.1, 2.1))
hist(dataframe_all$QUAL , breaks=100, col="cadetblue", border=F, main="" , 
     xlab="Qualité de variant (QUAL)", ylab = "Fréquence", xlim=c(0,max(dataframe_all$QUAL)))
dev.off()

# Distribution de QUAL entre b1 et b2
test_ <- wilcox.test(data_batch1$QUAL, data_batch2$QUAL, paired = FALSE)
p <- "< 2.2e-16"
png(file = "./2-images/18-Histogramme_et_boxplot_de_QUAL_entre_batch.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 5.5, 1.1, 2.1))
boxplot(dataframe_all$QUAL~dataframe_all$Batch, horizontal=TRUE, xaxt="n",
        col=c("brown1", "#00AFBB"), frame=F, cex = 0.3, pch=19, axes=FALSE, ann = FALSE)
par(mar=c(4, 5.5, 0, 2.1))
hist(as.numeric(data_batch2$QUAL), breaks=20, xlim=c(0,max(dataframe_all$QUAL)),
     col="#00AFBB", xlab="Qualité de variant (QUAL)", border=F,
     ylab="Fréquence", main="", cex.main = 2,
     cex.lab=3, ylim = c(0, 500000), cex.axis=2.5 )
hist(data_batch1$QUAL, breaks=20, xlim=c(0,max(dataframe_all$QUAL)), border=F,
     col="brown1", add=T)
# legend("right", legend=c("batch2","batch1"), col=c("#00AFBB", "brown1"),
#        pt.cex=4, pch=15, border = "white", bty = "n", cex = 3)
dev.off()

library(tidyverse)
library(viridis)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()
png(file = "./2-images/19-Violin_Chart_QUAL_entre_batch.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Batch, y=QUAL, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=30),
    axis.title.x = element_text(size=35, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=35, face="bold", vjust = 0.5),
    legend.title=element_blank(),
    legend.text=element_text(size=35),
    axis.text.x=element_text(size=35),
    axis.text.y=element_text(size=35)
  ) +
  xlab("Batch") + ylab("Qualité de variant (QUAL)")
dev.off()

# Distribution des QUAL par patient
png(file = "./2-images/20-Violin_Chart_QUAL_par_patient.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Patient, y=QUAL, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Patient") + ylab("Qualité de variant (QUAL)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()

# Distribution des QUAL par échantillon
png(file = "./2-images/21-Violin_Chart_QUAL_par_echantillon.png", width = 900, height = 800)
dataframe_all %>%
  ggplot( aes(x=Sample, y=QUAL, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(angle = 45)
  ) +
  xlab("Echantillon") + ylab("Qualité de variant (QUAL)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
#scale_color_viridis(discrete = TRUE, alpha=1, option="A")
#scale_color_brewer(palette="Dark2")
# scale_color_manual(values=c("#000000", "#FF0033", "#99CC00", "#FFFF00", "#996600",
#                             "#FF00CC", "#FF9933", "#0099FF", "#999999", "#990000",
#                             "#336600", "#CCCC33", "#663300", "#9900CC", "#0033FF"))
dev.off()
####################################################################################################################
# Distribution des mutations par Batch
mutation_per_batch <- table( dataframe_all$Batch )
mutation_per_batch_decreasing = sort( mutation_per_batch )
mutation_per_batch_decreasing_percentage <- prop.table( mutation_per_batch_decreasing )

mutation_per_batch_decreasing_percentage_df <- data.frame(mutation_per_batch_decreasing_percentage)
mutation_per_batch_decreasing_df <- data.frame(mutation_per_batch_decreasing)
df2 <- cbind(mutation_per_batch_decreasing_percentage_df, mutation_per_batch_decreasing_df$Freq)
colnames(df2) <- c("Batch", "Pourcentage", "Fréquence")
png(file = "./2-images/22-Avant_filtre_par_batch_ggplot.png", width = 900, height = 800)
ggplot(data=df2, aes(x=Batch, y=Pourcentage, fill=Batch)) +
  geom_bar(stat="identity")+
  theme_ipsum() +
  geom_text(aes(label=Fréquence), vjust=-1, color="Black", size=3.5)+
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    axis.text.x=element_text(size=35),
    axis.text.y=element_text(size=25),
    #legend.position = "none"
    legend.title=element_blank(),
    legend.text=element_text(size=35)
  ) +
  xlab("Batch") + ylab("Fréquence") + coord_flip()
dev.off()


png(file = "./2-images/22-Avant_filtre_par_batch_ggplot2.png", width = 900, height = 800)
library(scales)
ggplot(data=df2, aes(x="", y=Pourcentage, fill=Batch)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_ipsum() +
  geom_text(aes(y = Pourcentage/2 + c(0, cumsum(Pourcentage)[-length(Pourcentage)]), 
                label = percent(Pourcentage)), size=25) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title=element_blank(),
    #legend.text=element_text(size=35)
    axis.text.x=element_text(size=35),
    axis.text.y=element_text(size=25),
    legend.text=element_text(size=35),
    legend.position = "none"
  )
dev.off()
####################################################################################################################
# Distribution des mutations par patient
library(dplyr)
d1 <- distinct(dataframe_all, CHROM, POS, REF, ALT, .keep_all= TRUE)
mutation_per_patient <- table( d1$Patient )
mutation_per_patient_decreasing = sort( mutation_per_patient, decreasing=TRUE )
mutation_per_patient_decreasing_percentage <- prop.table( mutation_per_patient_decreasing )

mutation_per_patient_decreasing_percentage_df <- data.frame(mutation_per_patient_decreasing_percentage)
mutation_per_patient_decreasing_df <- data.frame(mutation_per_patient_decreasing)
df1 <- cbind(mutation_per_patient_decreasing_percentage_df, mutation_per_patient_decreasing_df$Freq)
colnames(df1) <- c("Patient", "Pourcentage", "Fréquence")

df1["Batch"] <- ifelse(as.numeric(levels(df1$Patient)) <= 6, "batch1", "batch2")
png(file = "./2-images/23-Avant_Filtre_par_patient_ggplot.png", width = 900, height = 800)
ggplot(data=df1, aes(x=Patient, y=Pourcentage, fill=Batch)) +
  geom_bar(stat="identity")+
  theme_ipsum() +
  #geom_text(aes(label=Fréquence), vjust=0, hjust=0.1, color="Black", size=6, angle=25) +
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    axis.text.x=element_text(size=35),
    axis.text.y=element_text(size=35),
    legend.title=element_blank(),
    legend.text=element_text(size=35)
  ) +
  ylab("Pourcentage")
dev.off()
####################################################################################################################
# Distribution des mutations par échantillon
mutation_per_sample <- table(dataframe_all$Sample)
mutation_per_sample_decreasing = sort( mutation_per_sample, decreasing = TRUE )
mutation_per_sample_decreasing_percentage <- prop.table( mutation_per_sample_decreasing )

mutation_per_sample_decreasing_percentage_df <- data.frame(mutation_per_sample_decreasing_percentage)
mutation_per_sample_decreasing_df <- data.frame(mutation_per_sample_decreasing)
df3 <- cbind(mutation_per_sample_decreasing_percentage_df, mutation_per_sample_decreasing_df$Freq)
colnames(df3) <- c("Sample", "Pourcentage", "Fréquence")

df3["Batch"] <- ifelse(substr(df3$Sample, 1, 3) == "D21", "batch2", "batch1")
png(file = "./2-images/24-Avant_filtre_par_echantillon.png", width = 900, height = 800)
ggplot(data=df3, aes(x=Sample, y=Pourcentage, fill=Batch)) +
  geom_bar(stat="identity") +
  theme_ipsum() +
  #geom_text(aes(label=Fréquence), hjust=0, color="Black", size=4) +
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    #axis.text.x = element_text(),
    axis.text.x=element_text(size=35, colour="black"),
    axis.text.y=element_text(size=23, colour="black"),
    #legend.position = "none"
    legend.title=element_blank(),
    legend.text=element_text(size=35)
  ) +
  xlab("Echantillon") + ylab("Pourcentage") + coord_flip()
dev.off()

#

mutation_per_sample <- table(dataframe_all$Sample)
mutation_per_sample_decreasing = sort( mutation_per_sample, decreasing = TRUE )
mutation_per_sample_decreasing_percentage <- prop.table( mutation_per_sample_decreasing )

mutation_per_sample_decreasing_percentage_df <- data.frame(mutation_per_sample_decreasing_percentage)
mutation_per_sample_decreasing_df <- data.frame(mutation_per_sample_decreasing)
df3 <- cbind(mutation_per_sample_decreasing_percentage_df, mutation_per_sample_decreasing_df$Freq)
colnames(df3) <- c("Sample", "Pourcentage", "Fréquence")

df3["Batch"] <- ifelse(substr(df3$Sample, 1, 3) == "D21", "batch2", "batch1")
png(file = "./2-images/24-Avant_filtre_par_echantillon_freq.png", width = 900, height = 800)
ggplot(data=df3, aes(x=Sample, y=Fréquence, fill=Batch)) +
  geom_bar(stat="identity") +
  theme_ipsum() +
  #geom_text(aes(label=Fréquence), hjust=0, color="Black", size=4) +
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    #axis.text.x = element_text(),
    axis.text.x=element_text(size=35, colour="black"),
    axis.text.y=element_text(size=23, colour="black"),
    #legend.position = "none"
    legend.title=element_blank(),
    legend.text=element_text(size=35)
  ) +
  #ylim(0, 600000) +
  scale_y_continuous(labels = comma_format(big.mark = "", #limits=c(0,600000),
                                           decimal.mark = ",")) +
  xlab("Echantillon") + ylab("Fréquence") + coord_flip()
dev.off()
####################################################################################################################
# Distribution des mutations par chromosome
print( "On calcule la distribution des mutations par chromosome" )
mutation_per_chr <- table( dataframe_all$CHROM )
mutation_per_chr_decreasing = sort( mutation_per_chr, decreasing = TRUE )
mutation_per_chr_decreasing_percentage <- prop.table( mutation_per_chr_decreasing )

mutation_per_chr_decreasing_percentage_df <- data.frame(mutation_per_chr_decreasing_percentage)
mutation_per_chr_decreasing_df <- data.frame(mutation_per_chr_decreasing)
df4 <- cbind(mutation_per_chr_decreasing_percentage_df, mutation_per_chr_decreasing_df$Freq)
colnames(df4) <- c("Chr", "Pourcentage", "Fréquence")

png(file = "./2-images/25-Avant_filtre_par_chr.png", width = 900, height = 800)
ggplot(data=df4, aes(x=Chr, y=Fréquence)) +
  geom_bar(stat="identity", fill="Cadetblue") +
  geom_text(aes(label=Fréquence), vjust=-1, color="Black", size=3) +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(angle = 45)
  ) +
  xlab("Chromosome") + ylab("Fréquence")
dev.off()
####################################################################################################################

#
write.table(dataframe_all, file="./1-text/Avant_filtre_avec_toutes_les_colonnes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

# Distribution des mutations après filtre
###########################################

# Importation des jeux de données
dataframe_all_patient_filtered <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/2-Filtre/2-Fichiers_filtrEs_avec_toutes_les_colonnes/echantillon_fusionnE-filtrE.txt", 
                                              header = FALSE, stringsAsFactors = FALSE  )

# On définit les noms des colonnes
colnames( dataframe_all_patient_filtered ) <- c( "CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "Details", "AF", "Sample", "Patient" )

# On crée une colonne DP
dataframe_all_patient_filtered["DP"] <- as.numeric( unlist( lapply( 1:nrow(dataframe_all_patient_filtered),
                                                                    function(i) { unlist( strsplit( dataframe_all_patient_filtered$Details[i], ":" ) )[3] } ) ) )

# On crée une colonne batch
dataframe_all_patient_filtered["Batch"] <- ifelse( (substr(dataframe_all_patient_filtered$Sample, 1, 3) == "D21"), "batch2", "batch1" )

# On crée une colonne AO
dataframe_all_patient_filtered["AO"] <- unlist(lapply( 1:nrow(dataframe_all_patient_filtered),
                                                       function(i) { as.numeric(unlist(strsplit(dataframe_all_patient_filtered$Details[i], ":"))[5]) } ))
####################################################################################################################
# On crée des un dataframe pour batch1 et un autre pour batch2
dataframe_all_patient_filtered$Patient <- as.character(dataframe_all_patient_filtered$Patient )
data_filtered_batch1 <- subset(dataframe_all_patient_filtered, Batch=="batch1")
data_filtered_batch2 <- subset(dataframe_all_patient_filtered, Batch=="batch2")

####################################################################################################################
# Distribution des DP
png(file = "./2-images/26-Histogramme_et_boxplot_de_DP.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 4, 1.1, 2.1))
boxplot(dataframe_all_patient_filtered$DP, horizontal=TRUE, ylim=c(0,max(dataframe_all_patient_filtered$DP)), xaxt="n",
        col=rgb(0.8,0.8,0,0.5), frame=F, cex = 0.3, pch=19)
par(mar=c(4, 4, 1.1, 2.1))
hist(dataframe_all_patient_filtered$DP , breaks=100, col="cadetblue", border=F, main="" , 
     xlab="Profondeur de lecture (DP)", ylab = "Fréquence", xlim=c(0,max(dataframe_all_patient_filtered$DP)))
dev.off()

# Distribution de DP entre b1 et b2
test_ <- wilcox.test(data_filtered_batch1$DP, data_filtered_batch2$DP, paired = FALSE)
p <- "= 2.328e-12"
png(file = "./2-images/27-Histogramme_et_boxplot_de_DP_entre_batch.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 5.5, 1.1, 2.1))
boxplot(log10(dataframe_all_patient_filtered$DP)~dataframe_all_patient_filtered$Batch, horizontal=TRUE, xaxt="n",
        col=c("brown1", "#00AFBB"), frame=F, cex = 0.3, pch=19, axes=FALSE, ann = FALSE,
        ylim=c(0,6))
par(mar=c(4, 5.5, 0, 2.1))
hist(log10(data_filtered_batch2$DP), breaks=30, xlim=c(0,6),
     col="#00AFBB", xlab="Profondeur de lecture (DP) en log10", border=F,
     ylab="Fréquence", main="",
     cex.main=2, cex.lab = 3, cex.axis=2.5, ylim = c(0, 5000) )
hist(log10(data_filtered_batch1$DP), breaks=20, xlim=c(0,6), border=F,
     col="brown1", add=T)
# legend("right", legend=c("Batch2","Batch1"), col=c("#00AFBB", "brown1"),
#        pt.cex=2, pch=15, border = "white", cex = 2, bty = "n")
dev.off()

png(file = "./2-images/28-Violin_Chart_DP_entre_batch.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Batch, y=DP, fill=Batch)) +
  geom_violin() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Profondeur de lecture (DP)")
dev.off()

# Distribution des DP par patient
png(file = "./2-images/29-Violin_Chart_DP_par_patient.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Patient, y=DP, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Patient") + ylab("Profondeur de lecture (DP)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()

# Distribution des DP par échantillon
png(file = "./2-images/30-Violin_Chart_DP_par_echantillon.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Sample, y=DP, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(angle = 45)
  ) +
  xlab("Echantillon") + ylab("Profondeur de lecture (DP)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()
####################################################################################################################
# Distribution des AF
png(file = "./2-images/31-Histogramme_et_boxplot_de_AF.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 4, 1.1, 2.1))
boxplot(dataframe_all_patient_filtered$AF, horizontal=TRUE, ylim=c(0,max(dataframe_all_patient_filtered$AF)), xaxt="n",
        col=rgb(0.8,0.8,0,0.5), frame=F, cex = 0.3, pch=19)
par(mar=c(4, 4, 1.1, 2.1))
hist(dataframe_all_patient_filtered$AF, breaks=100, col="cadetblue", border=F, main="" , 
     xlab="Fréquence allélique (AF)", ylab = "Fréquence", xlim=c(0,max(dataframe_all_patient_filtered$AF)))
dev.off()

# Distribution de AF entre b1 et b2
test_ <- wilcox.test(data_filtered_batch1$AF, data_filtered_batch2$AF, paired = FALSE)
p <- "< 2.2e-16"
png(file = "./2-images/32-Histogramme_et_boxplot_de_AF_entre_batch.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 5.5, 1.1, 2.1))
boxplot(dataframe_all_patient_filtered$AF~dataframe_all_patient_filtered$Batch, horizontal=TRUE, xaxt="n",
        col=c("brown1", "#00AFBB"), frame=F, cex = 0.3, pch=19, axes=FALSE, ann = FALSE)
par(mar=c(4, 5.5, 0, 2.1))
hist(data_filtered_batch2$AF, breaks=30, xlim=c(0,max(dataframe_all_patient_filtered$AF)),
     col="#00AFBB", xlab="Fréquence allélique (AF)", border=F,
     ylab="Fréquence", main="", cex.main = 2, cex.lab = 3, cex.axis = 2.5,
     ylim = c(0, 8000) )
hist(data_filtered_batch1$AF, breaks=30, xlim=c(0,max(dataframe_all_patient_filtered$AF)), border=F,
     col="brown1", add=T, cex.lab = 1.5)
legend("right", legend=c("batch2","batch1"), col=c("#00AFBB", "brown1"),
       pt.cex=5, pch=15, border = "white", bty = "n", cex = 3)
dev.off()

png(file = "./2-images/33-Violin_Chart_AF_entre_batch.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Batch, y=AF, fill=Batch)) +
  geom_violin() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Fréquence allélique (AF)")
dev.off()

# Distribution des AF par patient
png(file = "./2-images/34-Violin_Chart_AF_par_patient.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Patient, y=AF, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Patient") + ylab("Fréquence allélique (AF)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()

# Distribution des AF par échantillon
png(file = "./2-images/35-Violin_Chart_AF_par_echantillon.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Sample, y=AF, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(angle = 45)
  ) +
  xlab("Echantillon") + ylab("Fréquence allélique (AF)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()
####################################################################################################################
# Distribution des AF selon DP par batch
png(file =  "./2-images/36-Nuage_de_points_de_AF_DP_pour_batch1.png", width = 900, height = 800)
plot( data_filtered_batch1$DP, data_filtered_batch1$AF, main = "Batch 1",
      pch=19, xlab = "Profondeur de lecture (DP)", ylab = "Fréquence allélique (AF)", cex = 0.5)
dev.off()

png(file =  "./2-images/37-Nuage_de_points_de_AF_DP_pour_batch2.png", width = 900, height = 800)
plot( data_filtered_batch2$DP, data_filtered_batch2$AF, main = "Batch 2",
      pch=19, xlab = "Profondeur de lecture (DP)", ylab = "Fréquence allélique (AF)", cex = 0.5)
dev.off()

####

png(file =  "./2-images/38-Nuage_de_points_de_AF_DP_pour_batch1_log_y.png", width = 900, height = 800)
plot( data_filtered_batch1$DP, data_filtered_batch1$AF, main = "Batch 1", log = "y",
      pch=19, xlab = "Profondeur de lecture (DP)", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
dev.off()

png(file =  "./2-images/39-Nuage_de_points_de_AF_DP_pour_batch2_log_y.png", width = 900, height = 800)
plot( data_filtered_batch2$DP, data_filtered_batch2$AF, main = "Batch 2", log = "y",
      pch=19, xlab = "Profondeur de lecture (DP)", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
dev.off()

####

png(file =  "./2-images/40-Nuage_de_points_de_AF_DP_pour_batch1_log.png", width = 900, height = 800)
plot( data_filtered_batch1$DP, data_filtered_batch1$AF, main = "Batch 1", log = "xy",
      pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
dev.off()

png(file =  "./2-images/41-Nuage_de_points_de_AF_DP_pour_batch2_log.png", width = 900, height = 800)
plot( data_filtered_batch2$DP, data_filtered_batch2$AF, main = "Batch 2", log = "xy",
      pch=19, xlab = "Profondeur de lecture (log(DP))", ylab = "Fréquence allélique (log(AF))", cex = 0.5)
dev.off()
####################################################################################################################
# Distribution des QUAL
png(file = "./2-images/42-Histogramme_et_boxplot_de_QUAL.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 4, 1.1, 2.1))
boxplot(dataframe_all_patient_filtered$QUAL, horizontal=TRUE, ylim=c(0,max(dataframe_all_patient_filtered$QUAL)), xaxt="n",
        col=rgb(0.8,0.8,0,0.5), frame=F, cex = 0.3, pch=19)
par(mar=c(4, 4, 1.1, 2.1))
hist(dataframe_all_patient_filtered$QUAL , breaks=100, col="cadetblue", border=F, main="" , 
     xlab="Qualité de variant (QUAL)", ylab = "Fréquence", xlim=c(0,max(dataframe_all_patient_filtered$QUAL)))
dev.off()

# Distribution de QUAL entre b1 et b2
test_ <- wilcox.test(data_filtered_batch1$DP, data_filtered_batch2$DP, paired = FALSE)
p <- "= 2.328e-12"
png(file = "./2-images/43-Histogramme_et_boxplot_de_QUAL_entre_batch.png", width = 900, height = 800)
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 5.5, 1.1, 2.1))
boxplot(dataframe_all_patient_filtered$QUAL~dataframe_all_patient_filtered$Batch, horizontal=TRUE, xaxt="n",
        col=c("brown1", "#00AFBB"), frame=F, cex = 0.3, pch=19, axes=FALSE, ann = FALSE)
par(mar=c(4, 5.5, 0, 2.1))
hist(data_filtered_batch2$QUAL, breaks=20, xlim=c(0,max(dataframe_all_patient_filtered$QUAL)),
     col="#00AFBB", xlab="Qualité de variant (QUAL)", border=F,
     ylab="Fréquence", main="", cex.lab=3, cex.axis=2.5, 
     cex.main = 2, ylim = c(0,20000) )
hist(data_filtered_batch1$QUAL, breaks=20, xlim=c(0,max(dataframe_all_patient_filtered$QUAL)), border=F,
     col="brown1", add=T)
# legend("right", legend=c("Batch2","Batch1"), col=c("#00AFBB", "brown1"),
#        pt.cex=2, pch=15, border = "white", bty = "n", cex=2)
dev.off()

png(file = "./2-images/44-Violin_Chart_QUAL_entre_batch.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Batch, y=QUAL, fill=Batch)) +
  geom_violin() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Batch") + ylab("Qualité de variant (QUAL)")
dev.off()

# Distribution des QUAL par patient
png(file = "./2-images/45-Violin_Chart_QUAL_par_patient.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Patient, y=QUAL, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  xlab("Patient") + ylab("Qualité de variant (QUAL)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()

# Distribution des QUAL par échantillon
png(file = "./2-images/46-Violin_Chart_QUAL_par_echantillon.png", width = 900, height = 800)
dataframe_all_patient_filtered %>%
  ggplot( aes(x=Sample, y=QUAL, fill=Batch)) +
  geom_violin() +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(angle = 45)
  ) +
  xlab("Echantillon") + ylab("Qualité de variant (QUAL)") +
  scale_fill_manual(values=c("brown1", "#00AFBB"))
dev.off()
####################################################################################################################
# Distribution des mutations par Batch
mutation_per_batch <- table( dataframe_all_patient_filtered$Batch )
mutation_per_batch_decreasing = sort( mutation_per_batch )
mutation_per_batch_decreasing_percentage <- prop.table( mutation_per_batch_decreasing )

mutation_per_batch_decreasing_percentage_df <- data.frame(mutation_per_batch_decreasing_percentage)
mutation_per_batch_decreasing_df <- data.frame(mutation_per_batch_decreasing)
df2 <- cbind(mutation_per_batch_decreasing_percentage_df, mutation_per_batch_decreasing_df$Freq)
colnames(df2) <- c("Batch", "Pourcentage", "Fréquence")
png(file = "./2-images/47-Apres_filtre_par_batch_ggplot.png", width = 900, height = 800)
ggplot(data=df2, aes(x=Batch, y=Pourcentage, fill=c("brown1", "#00AFBB"))) +
  geom_bar(stat="identity")+
  theme_ipsum() +
  geom_text(aes(label=Fréquence), vjust=-1, color="Black", size=3.5)+
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    #axis.text.x = element_text(),
    axis.text.x=element_text(size=35),
    axis.text.y=element_text(size=25),
    #legend.position = "none"
    legend.title=element_blank(),
    legend.text=element_text(size=35)
  ) +
  xlab("Batch") + ylab("Pourcentage")
dev.off()

png(file = "./2-images/47-Apres_filtre_par_batch_ggplot2.png", width = 900, height = 800)
library(scales)
ggplot(data=df2, aes(x="", y=Pourcentage, fill=c("brown1", "#00AFBB"))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_ipsum() +
  geom_text(aes(y = Pourcentage/2 + c(0, cumsum(Pourcentage)[-length(Pourcentage)]),
                label = percent(Pourcentage)), size=25) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title=element_blank(),
    #legend.text=element_text(size=35)
    axis.text.x=element_text(size=35),
    axis.text.y=element_text(size=25),
    legend.text=element_text(size=35),
    legend.position = "none"
  )
dev.off()

png(file = "./2-images/47-Apres_filtre_par_batch_r.png", width = 900, height = 800)
pie(table(dataframe_all_patient_filtered$Batch), col=c("brown1", "#00AFBB"), labels="")
legend("topright", legend=c("Batch1","Batch2"), col=c("brown1", "#00AFBB"),
       pt.cex=2, pch=15, border = "white", bty = "n", cex=2)
dev.off()
####################################################################################################################
# Distribution des mutations par patient
d2 <- distinct(dataframe_all_patient_filtered, CHROM, POS, REF, ALT, .keep_all= TRUE)
mutation_per_patient <- table( d2$Patient )
mutation_per_patient_decreasing = sort( mutation_per_patient, decreasing=TRUE )
mutation_per_patient_decreasing_percentage <- prop.table( mutation_per_patient_decreasing )

mutation_per_patient_decreasing_percentage_df <- data.frame(mutation_per_patient_decreasing_percentage)
mutation_per_patient_decreasing_df <- data.frame(mutation_per_patient_decreasing)
df1 <- cbind(mutation_per_patient_decreasing_percentage_df, mutation_per_patient_decreasing_df$Freq)
colnames(df1) <- c("Patient", "Pourcentage", "Fréquence")

df1["Batch"] <- ifelse(as.numeric(levels(df1$Patient)) <= 6, "batch1", "batch2")
png(file = "./2-images/48-Apres_Filtre_par_patient_ggplot.png", width = 900, height = 800)
ggplot(data=df1, aes(x=Patient, y=Pourcentage, fill=Batch)) +
  geom_bar(stat="identity") +
  theme_ipsum() +
  geom_text(aes(label=Fréquence), vjust=-1, color="Black", size=3.5) +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
dev.off()
####################################################################################################################
# Distribution des mutations par échantillon
mutation_per_sample <- table(dataframe_all_patient_filtered$Sample)
mutation_per_sample_decreasing = sort( mutation_per_sample, decreasing = TRUE )
mutation_per_sample_decreasing_percentage <- prop.table( mutation_per_sample_decreasing )

mutation_per_sample_decreasing_percentage_df <- data.frame(mutation_per_sample_decreasing_percentage)
mutation_per_sample_decreasing_df <- data.frame(mutation_per_sample_decreasing)
df3 <- cbind(mutation_per_sample_decreasing_percentage_df, mutation_per_sample_decreasing_df$Freq)
colnames(df3) <- c("Sample", "Pourcentage", "Fréquence")

df3["Batch"] <- ifelse(substr(df3$Sample, 1, 3) == "D21", "batch2", "batch1")
png(file = "./2-images/49-Apres_filtre_par_echantillon.png", width = 900, height = 800)
ggplot(data=df3, aes(x=Sample, y=Pourcentage, fill=Batch)) +
  geom_bar(stat="identity") +
  theme_ipsum() +
  #geom_text(aes(label=Fréquence), hjust=-0.2, color="Black", size=3.5) +
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    #axis.text.x = element_text(),
    axis.text.x=element_text(size=35, colour="black"),
    axis.text.y=element_text(size=25, colour="black"),
    #legend.position = "none"
    legend.title=element_blank(),
    legend.text=element_text(size=35)
  ) +
  xlab("Echantillon") + ylab("Pourcentage") + coord_flip()
dev.off()

#

mutation_per_sample <- table(dataframe_all_patient_filtered$Sample)
mutation_per_sample_decreasing = sort( mutation_per_sample, decreasing = TRUE )
mutation_per_sample_decreasing_percentage <- prop.table( mutation_per_sample_decreasing )

mutation_per_sample_decreasing_percentage_df <- data.frame(mutation_per_sample_decreasing_percentage)
mutation_per_sample_decreasing_df <- data.frame(mutation_per_sample_decreasing)
df3 <- cbind(mutation_per_sample_decreasing_percentage_df, mutation_per_sample_decreasing_df$Freq)
colnames(df3) <- c("Sample", "Pourcentage", "Fréquence")

df3["Batch"] <- ifelse(substr(df3$Sample, 1, 3) == "D21", "batch2", "batch1")
png(file = "./2-images/49-Apres_filtre_par_echantillon.png", width = 900, height = 800)
ggplot(data=df3, aes(x=Sample, y=Fréquence, fill=Batch)) +
  geom_bar(stat="identity") +
  theme_ipsum() +
  #geom_text(aes(label=Fréquence), hjust=-0.2, color="Black", size=3.5) +
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    #axis.text.x = element_text(),
    axis.text.x=element_text(size=35, colour="black"),
    axis.text.y=element_text(size=25, colour="black"),
    #legend.position = "none"
    legend.title=element_blank(),
    legend.text=element_text(size=35)
  ) +
  xlab("Echantillon") + ylab("Fréquence") + coord_flip()
dev.off()
####################################################################################################################
# Distribution des mutations par chromosome
print( "On calcule la distribution des mutations par chromosome" )
mutation_per_chr <- table( dataframe_all_patient_filtered$CHROM )
mutation_per_chr_decreasing = sort( mutation_per_chr, decreasing = TRUE )
mutation_per_chr_decreasing_percentage <- prop.table( mutation_per_chr_decreasing )

mutation_per_chr_decreasing_percentage_df <- data.frame(mutation_per_chr_decreasing_percentage)
mutation_per_chr_decreasing_df <- data.frame(mutation_per_chr_decreasing)
df4 <- cbind(mutation_per_chr_decreasing_percentage_df, mutation_per_chr_decreasing_df$Freq)
colnames(df4) <- c("Chr", "Pourcentage", "Fréquence")

png(file = "./2-images/50-Apres_filtre_par_chr.png", width = 900, height = 800)
ggplot(data=df4, aes(x=Chr, y=Fréquence)) +
  geom_bar(stat="identity", fill="Cadetblue") +
  geom_text(aes(label=Fréquence), vjust=-1, color="Black", size=3) +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x = element_text(angle = 45)
  ) +
  xlab("Chromosome") + ylab("Fréquence")
dev.off()
####################################################################################################################

#
write.table(dataframe_all_patient_filtered, file="./1-text/Apres_filtre_avec_toutes_les_colonnes.txt", col.names=T, quote=F, row.names=F, sep = "\t")
####################################################################################################################

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

# Distribution des mutations après annotation
#############################################

# Importation des données
Annotated <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/3-Annotation/2-Fichiers_annotEs/echantillon_fusionnE-annotE.txt",
                         stringsAsFactors = FALSE, sep = "" )

colnames(Annotated) <- c("Uploaded_variation",	"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence",
                         "cDNA_position",	"CDS_position",	"Protein_position",	"Amino_acids",	"Codons",	"Existing_variation",	"Extra")

# Distribution des mutations par gène driver
print( "Importation du jeu de données GENES-DRIVER" )
drivers <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/0-Metadata/IntOGen-DriverGenes_HNSC_2022.tsv", header = TRUE, sep = "\t", dec = ".",
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
write.table(Annotated_reduced, file="./1-text/Annotation_reduit_avec_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

# On calcule le nombre de mutations pour chaque gène driver dans nos jeux de données
print( "On calcule le nombre de mutations pour chaque gène driver dans nos jeux de données" )
driver_list <- data.frame( Symbol = drivers$Symbol, All_mutations = rep( NA, length(drivers$Symbol) ), stringsAsFactors = FALSE )

for (i in 1:nrow(driver_list)) {
  subset_i <- subset( Annotated_reduced, Gene==drivers$Symbol[i] )
  if (nrow(subset_i)!=0) {
    driver_list$All_mutations[i] <- nrow(subset_i)
  }
  else
  {
    driver_list$All_mutations[i] <- 0
  }
}
driver_list_sorted <- driver_list[order(as.numeric(driver_list$All_mutations), decreasing = TRUE), ]
write.table(driver_list_sorted, file="./1-text/Nombre_de_mutations_par_gene_driver.txt", col.names=T, quote=F, row.names=F, sep = "\t")

# Les 500 gènes les plus mutés parmi les gènes exprimés dans les cellules épithéliales de la muqueuse orale
print( "On cherche les 500 gènes les plus mutés parmi les gènes exprimés dans les cellules épithéliales de la muqueuse orale" )

expressed_genes <- read.table( file = "../4-Matrice_atlas/1-Sortie/genes_exprimes_dans_les_cellulles_epitheliales_de_la_muqueuse_orale.txt",
                               stringsAsFactors = FALSE )
colnames(expressed_genes) <- "Genes"
expressed_genes["Mutations_number"] <- rep("NA", nrow(expressed_genes))
for (i in 1:nrow(expressed_genes)) {
  subset_i <- subset( Annotated_reduced, Gene==expressed_genes$Genes[i] )
  if (nrow(subset_i)!=0) {
    expressed_genes$Mutations_number[i] <- nrow(subset_i)
  }
  else
  {
    expressed_genes$Mutations_number[i] <- 0
  }
}
expressed_genes <- expressed_genes[!(expressed_genes$Genes %in% driver_list$Symbol),]
expressed_genes_sorted <- expressed_genes[order(as.numeric(expressed_genes$Mutations_number), decreasing = TRUE), ]
write.table(expressed_genes_sorted, file="./1-text/Genes_exprimes_et_nombre_de_mutations.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table( data.frame(expressed_genes_sorted[1:500,], stringsAsFactors = FALSE), 
             file="./1-text/500_genes_les_plus_mutes_parmi_les_genes_exprimes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

# Les 1000 gènes les plus mutés parmi les gènes non exprimés dans les cellules épithéliales de la muqueuse orale
print( "On cherche les 1000 gènes les plus mutés parmi les gènes non exprimés dans les cellules épithéliales de la muqueuse orale" )

unexpressed_genes <- read.table( file = "../4-Matrice_atlas/1-Sortie/genes_non_exprimes_dans_les_cellulles_epitheliales_de_la_muqueuse_orale.txt",
                                 stringsAsFactors = FALSE )
colnames(unexpressed_genes) <- "Genes"
unexpressed_genes["Mutations_number"] <- rep("NA", nrow(unexpressed_genes))
for (i in 1:nrow(unexpressed_genes)) {
  subset_i <- subset( Annotated_reduced, Gene==unexpressed_genes$Genes[i] )
  if (nrow(subset_i)!=0) {
    unexpressed_genes$Mutations_number[i] <- nrow(subset_i)
  }
  else
  {
    unexpressed_genes$Mutations_number[i] <- 0
  }
}
unexpressed_genes <- unexpressed_genes[!(unexpressed_genes$Genes %in% driver_list$Symbol),]
unexpressed_genes_sorted <- unexpressed_genes[order(as.numeric(unexpressed_genes$Mutations_number), decreasing = TRUE), ]
write.table(unexpressed_genes_sorted, file="./1-text/Genes_non_exprimes_et_nombre_de_mutations.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table( data.frame(unexpressed_genes_sorted[1:1000,], stringsAsFactors = FALSE), 
             file="./1-text/1000_genes_les_plus_mutes_parmi_les_genes_non_exprimes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

#############################################################################################################################################################

