rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes/Analyse")

ids <- c("D181208", "D181210", "D181211", "D181212", "D181213", "D181215", "D202061", "D202062",
         "D202063", "D202064", "D202065", "D202066", "D210285", "D210288", "D210289", "D210290",
         "D210294", "D210295", "D210296", "D210297", "D210327", "D210328", "D210332", "D210334",
         "D210335", "D210338", "D210339")

annotation <- read.table("~/ZINARA/Gnomic/1-StageM2/5-Stats/1-text/Annotation_reduit_avec_genes.txt",
                         header = TRUE, sep = "\t")
annotation <- annotation[,1:14]
annotation["Uploaded_variation"] <- unlist(lapply(1:nrow(annotation), function(i) {paste(unlist(strsplit(annotation$Uploaded_variation[i], "_"))[1],
                                                                 unlist(strsplit(annotation$Uploaded_variation[i], "_"))[2],
                                                                 unlist(strsplit(annotation$Uploaded_variation[i], "_"))[3],
                                                                 unlist(strsplit(annotation$Uploaded_variation[i], "_"))[4],
                                                                 unlist(strsplit(annotation$Uploaded_variation[i], "_"))[5],
                                                                 sep = "_")}))

expressed <- read.table("~/ZINARA/Gnomic/1-StageM2/4-Matrice_atlas/1-Sortie/genes_exprimes_dans_les_cellulles_epitheliales_de_la_muqueuse_orale.txt",
                        header = FALSE, sep = "\t")

neo_Ag_summary <- data.frame("Sample"=character(),"Total"=integer(),"Total_WB"=integer(),"Total_SB"=integer(),
                             "Total_Region_0"=integer(),"Total_WB_Region_0"=integer(),"Total_SB_Region_0"=integer(),
                             "Clonal"=integer(),"Subclonal"=integer(),"Shared"=integer(),"Clonal_WB"=integer(),
                             "Clonal_SB"=integer(),"Subclonal_WB"=integer(),"Subclonal_SB"=integer(),"Shared_WB"=integer(),
                             "Shared_SB"=integer())
for (i in 1:2) {
  neo_Ag_summary_i <- read.table(sprintf("~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes/NeoPredPipe/Sortie/b1_%s/TestRun.neoantigens.summarytable.txt", i), 
                                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  neo_Ag_summary <- rbind(neo_Ag_summary, neo_Ag_summary_i)
}
for (i in 1:2) {
  neo_Ag_summary_i <- read.table(sprintf("~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes/NeoPredPipe/Sortie/b2_%s/TestRun.neoantigens.summarytable.txt", i), 
                                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  neo_Ag_summary <- rbind(neo_Ag_summary, neo_Ag_summary_i)
}

neo_Ag <- data.frame("Sample"=character(),"Region"=integer(),"Line"=character(),"chr"=character(),"allelepos"=integer(),"ref"=character(),
                     "alt"=character(),"GeneName:RefSeqID"=character(),"pos"=integer(),"hla"=character(),"peptide"=character(),
                     "core"=character(),"Of"=integer(),"Gp"=integer(),"Gl"=integer(),"Ip"=integer(),"Il"=integer(),"Icore"=character(),
                     "Identity"=character(),"Score"=numeric(),"Binding Affinity"=numeric(),"Rank"=numeric(),"Candidate"=character(),
                     "BindLevel"=character(),"Novelty"=integer(),"id"=character())
for (i in 1:2) {
  neo_Ag_i <- read.table(sprintf("~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes/NeoPredPipe/Sortie/b1_%s/TestRun.neoantigens.txt", i), 
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(neo_Ag_i) <- c("Sample","Region","Line","chr","allelepos","ref","alt","GeneName.RefSeqID","pos","hla","peptide","core",
                        "Of","Gp","Gl","Ip","Il","Icore","Identity","Score","Binding.Affinity","Rank","Candidate","BindLevel","Novelty")
  neo_Ag <- rbind(neo_Ag, neo_Ag_i)
}
for (i in 1:2) {
  neo_Ag_i <- read.table(sprintf("~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes/NeoPredPipe/Sortie/b2_%s/TestRun.neoantigens.txt", i), 
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(neo_Ag_i) <- c("Sample","Region","Line","chr","allelepos","ref","alt","GeneName.RefSeqID","pos","hla","peptide","core",
                          "Of","Gp","Gl","Ip","Il","Icore","Identity","Score","Binding.Affinity","Rank","Candidate","BindLevel","Novelty")
  neo_Ag <- rbind(neo_Ag, neo_Ag_i)
}
neo_Ag["id"] <- paste(neo_Ag$chr, neo_Ag$allelepos, neo_Ag$ref, neo_Ag$alt, neo_Ag$Sample, sep = "_")
neo_Ag["Gene"] <- unlist(lapply(1:nrow(neo_Ag), function(i) {unlist(strsplit(neo_Ag$GeneName.RefSeqID[i], ":"))[1]}))

# Seulement les mutations sur les gènes exprimés dans l'épithélium de la muqueuse orale
neo_Ag_expressed <- neo_Ag[neo_Ag$Gene %in% expressed$V1,]

df_all <- data.frame(Sample=character(), Mutation=numeric(), MutationGivingNeoantigenWithinExpressedGenes=numeric(), 
                     p_MutationGivingNeoantigenWithinExpressedGenes=numeric(), NeoantigenTotal=numeric(), NeoantigenByExpressedGenes=numeric(),
                     NewNeoantigenByExpressedGenes=numeric(), p_NewNeoantigenByExpressedGenes=numeric(), NewNeoantigenWBByExpressedGenes=numeric(),
                     p_NewNeoantigenWBByExpressedGenes_among_NewNeoantigen=numeric(), p_NewNeoantigenWBByExpressedGenes_among_Neoantigen=numeric(), NewNeoantigenSBByExpressedGenes=numeric(),
                     p_NewNeoantigenSBByExpressedGenes_among_NewNeoantigen=numeric(), p_NewNeoantigenSBByExpressedGenes_among_Neoantigen=numeric())

for (idx in ids) {
  
  data_i <- read.table(sprintf("~/ZINARA/Gnomic/1-StageM2/2-Filtre/3-Fichiers_filtrEs_vcf/%s.vcf", idx), 
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(data_i) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
  data_i["ID"] <- unlist(lapply(1:nrow(data_i), function(i) {paste(unlist(strsplit(data_i$ID[i], "_"))[1],
                                                                   unlist(strsplit(data_i$ID[i], "_"))[2],
                                                                   unlist(strsplit(data_i$ID[i], "_"))[3],
                                                                   unlist(strsplit(data_i$ID[i], "_"))[4],
                                                                   unlist(strsplit(data_i$ID[i], "_"))[5],
                                                                   sep = "_")}))
  
  neo_Ag_expressed_i <- subset(neo_Ag_expressed, Sample==idx)
  neo_Ag_summary_i <- subset(neo_Ag_summary, Sample==idx)
  
  id_data_i <- unique(data_i$ID)
  id_neo_Ag_expressed_i <- unique(neo_Ag_expressed_i$id)
  intersection_i <- intersect(id_data_i, id_neo_Ag_expressed_i)
  difference_i <- setdiff(id_data_i, id_neo_Ag_expressed_i)
  
  # Pourcentage de mutations donnant lieu à des néoantigènes
  p_i <- round(length(intersection_i) / length(id_data_i) * 100, 2)
  # Nombre de néoantigènes exprimés
  N <- nrow(neo_Ag_expressed_i)
  # Nombre total de néoantigènes
  N_ <- nrow(subset(neo_Ag, Sample==idx))
  # Pourcentage de nouveaux néoantigènes exprimés
  new <- nrow(subset(neo_Ag_expressed_i, Novelty==1))
  p_new <- round( new / N_ * 100, 2)
  # Pourcentage des nouveaux néoantigènes dans les gènes exprimés qui sont des ligants fort parmi les nouveaux néoantigènes
  SB_i <- nrow(subset(neo_Ag_expressed_i, Novelty==1 & BindLevel=="SB"))
  posx <- grep(idx, neo_Ag_summary_i$Sample)
  p_SB_i <- round( SB_i / neo_Ag_summary_i$Total[posx] * 100, 2)
  # Pourcentage des nouveaux néoantigènes qui sont des ligants fort parmi tous les néoantigènes
  p_SB_i_ <- round( SB_i / N_ * 100, 2)
  # Pourcentage des nouveaux néoantigènes dans les gènes exprimés qui sont des ligants faible parmi les nouveaux néoantigènes
  WB_i <- nrow(subset(neo_Ag_expressed_i, Novelty==1 & BindLevel=="WB"))
  p_WB_i <- round( WB_i / neo_Ag_summary_i$Total[posx] * 100, 2)
  # Pourcentage des nouveaux néoantigènes qui sont des ligants faible parmi tous les néoantigènes
  p_WB_i_ <- round( WB_i / N_ * 100, 2)
  
  df_all[dim(df_all)[1]+1,] <- c(idx, length(id_data_i), length(intersection_i), p_i, N_, N, new, p_new, WB_i,
                                 p_WB_i, p_WB_i_, SB_i, p_SB_i, p_SB_i_)
}
write.table(neo_Ag_expressed, file="./sortie_exprimE/Neoantigene_exprimE.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table(neo_Ag_summary, file="./sortie_exprimE/Neoantigene_synthese.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table(df_all, file="./sortie_exprimE/Resultat.txt", col.names=T, quote=F, row.names=F, sep = "\t")


### PLOT 1 : Mutations donnant lieu à des néoantigènes parmi toutes les mutations par échantillon
data_2 <- read.table("~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes/Analyse/sortie/Resultat.txt", header = TRUE,
                     sep = "\t", stringsAsFactors = FALSE)
d1 <- df_all[, c("Sample", "p_MutationGivingNeoantigenWithinExpressedGenes")]
d1$p_MutationGivingNeoantigenWithinExpressedGenes <- as.numeric(d1$p_MutationGivingNeoantigenWithinExpressedGenes)
d1["p_MutationGivingNeoantigenWithinUnexpressedGenes"] <- data_2$p_MutationGivingNeoantigen - d1$p_MutationGivingNeoantigenWithinExpressedGenes
d1["NoNeoAg"] <- 100-d1$p_MutationGivingNeoantigenWithinExpressedGenes-d1$p_MutationGivingNeoantigenWithinUnexpressedGenes
td1 <- t(d1)
td1 <- data.frame(td1)
colnames(td1) <- d1$Sample
td1 <- td1[2:4,]
td1 <- as.matrix(td1)
png(file = "./sortie_exprimE/images/1-Mutations_donnant_neoantigenes.png", width = 900, height = 1100)
par(mar=c(5, 6, 5, 1))
barplot(td1, las=1, col=c("#5F9EA0", "#E1B378", "grey"), xlab = "Pourcentage", cex.lab=3, cex.axis=3,
        horiz = TRUE)
legend("right", legend=c("Mutations donnant des neoAg \n chez les gènes exprimés \n","Mutations donnant des neoAg \n chez les gènes non exprimés","Mutations ne donnant pas de neoAg"), col=c("#5F9EA0", "#E1B378", "grey"),
       pt.cex=2, pch=15, border = "white", cex = 2.5, bty = "o")
dev.off()

#####################
n = 3
Echantillon <- c(rep("D181208", n), rep("D181210", n), rep("D181211", n), rep("D181212", n),
                 rep("D181213", n), rep("D181215", n), rep("D202061", n), rep("D202062", n),
                 rep("D202063", n), rep("D202064", n), rep("D202065", n), rep("D202066", n),
                 rep("D210285", n), rep("D210288", n), rep("D210289", n), rep("D210290", n),
                 rep("D210294", n), rep("D210295", n), rep("D210296", n), rep("D210297", n),
                 rep("D210327", n), rep("D210328", n), rep("D210332", n), rep("D210334", n),
                 rep("D210335", n), rep("D210338", n), rep("D210339", n)
)
Categorie <- rep(c("Mutations donnant \n des neoAg \n chez les gènes \n exprimés",
                   "Mutations donnant \n des neoAg \n chez les gènes \n non exprimés",
                   "Mutations \n ne donnant pas \n de neoAg"), 27)
Pourcentage <- c()
for (i in 1:27) {
  for (j in 1:n) {
    Pourcentage <- c(Pourcentage, as.numeric(td1[j,i]))
  }
}
d1 <- data.frame(Echantillon, Categorie, Pourcentage)

library(tidyverse)
library(viridis)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()

png(file = "./sortie_exprimE/images/1-Mutations_donnant_neoantigenes_ggplot.png", width = 900, height = 1000)
ggplot(d1, aes(fill = Categorie, y = Pourcentage, x = Echantillon))+
  theme_ipsum()+
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("")+
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    #axis.text.x = element_text(),
    axis.text.x=element_text(size=35),
    axis.text.y=element_text(size=25),
    #legend.position = "none"
    legend.title=element_blank(),
    legend.text=element_text(size=30),
    legend.position="top"
  )+ coord_flip()
dev.off()
#####################

### PLOT 2 : Nouveaux néoantigènes parmi tout les néoantigènes par échantillon
d2 <- df_all[, c("Sample", "p_NewNeoantigenWBByExpressedGenes_among_Neoantigen", "p_NewNeoantigenSBByExpressedGenes_among_Neoantigen")]
d2["p_NewNeoantigenByNunexpressedGenes_among_Neoantigen"] <- round((as.numeric(data_2$NewNeoantigen )-as.numeric(df_all$NewNeoantigenByExpressedGenes))/as.numeric(data_2$Neoantigen) * 100, 2)
d2["p_knownNeoantigene"] <- round((as.numeric(data_2$Neoantigen)-as.numeric(data_2$NewNeoantigen))/as.numeric(data_2$Neoantigen)*100, 2)
td2 <- t(d2)
td2 <- data.frame(td2)
colnames(td2) <- d2$Sample
td2 <- td2[2:5,]
td2 <- as.matrix(td2)
png(file = "./sortie_exprimE/images/2-Nouveaux_neoantigenes_parmi_neoantigenes_total.png", width = 900,
    height = 1100)
par(mar=c(5, 6, 5, 1))
barplot(td2, las=1, col=c("#5F9EA0", "#E1B378", "#D6604D", "#D1E5F0"), xlab = "Pourcentage", cex.lab=3, cex.axis=3,
        horiz = TRUE)
legend("center", legend=c("Nouveaux neoAg 'ligant faible' \n chez les gènes exprimés \n","Nouveaux neoAg 'ligant fort' \n chez les gènes exprimés \n",
                         "Nouveaux neoAg chez les gènes \n non exprimés", "neoAg connus"), col=c("#5F9EA0", "#E1B378", "#D6604D", "#D1E5F0"),
       pt.cex=2, pch=15, border = "white", cex = 2.5, bty = "o")
dev.off()

#####################
n = 4
Echantillon <- c(rep("D181208", n), rep("D181210", n), rep("D181211", n), rep("D181212", n),
                 rep("D181213", n), rep("D181215", n), rep("D202061", n), rep("D202062", n),
                 rep("D202063", n), rep("D202064", n), rep("D202065", n), rep("D202066", n),
                 rep("D210285", n), rep("D210288", n), rep("D210289", n), rep("D210290", n),
                 rep("D210294", n), rep("D210295", n), rep("D210296", n), rep("D210297", n),
                 rep("D210327", n), rep("D210328", n), rep("D210332", n), rep("D210334", n),
                 rep("D210335", n), rep("D210338", n), rep("D210339", n)
)
Categorie <- rep(c("Nouveaux neoAg \n 'ligant faible' \n chez les gènes \n exprimés",
                   "Nouveaux neoAg \n 'ligant fort' \n chez les gènes \n exprimés",
                   "Nouveaux neoAg \n chez les gènes \n non exprimés", "neoAg \n connus"), 27)
Pourcentage <- c()
for (i in 1:27) {
  for (j in 1:n) {
    Pourcentage <- c(Pourcentage, as.numeric(td2[j,i]))
  }
}
d2 <- data.frame(Echantillon, Categorie, Pourcentage)
png(file = "./sortie_exprimE/images/2-Nouveaux_neoantigenes_parmi_neoantigenes_total_ggplot.png", width = 900, height = 1000)
ggplot(d2, aes(fill = Categorie, y = Pourcentage, x = Echantillon))+
  theme_ipsum()+
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("")+
  theme(
    plot.title = element_text(size=35),
    axis.title.x = element_text(size=35, face="bold"),
    axis.title.y = element_text(size=35, face="bold"),
    #axis.text.x = element_text(),
    axis.text.x=element_text(size=35),
    axis.text.y=element_text(size=25),
    #legend.position = "none"
    legend.title=element_blank(),
    legend.text=element_text(size=30),
    legend.position="top"
  )+ coord_flip()
dev.off()