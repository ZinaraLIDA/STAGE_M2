rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/10-Evolution")
#############################################################
patients <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/0-Metadata/patients.txt", header = TRUE,
                       sep = "\t")
#############################################################
Annotated <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/3-Annotation/2-Fichiers_annotEs/echantillon_fusionnE-annotE.txt",
                        header = FALSE, sep = "\t")
colnames(Annotated) <- c("Uploaded_variation",	"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence",
                         "cDNA_position",	"CDS_position",	"Protein_position",	"Amino_acids",	"Codons",	"Existing_variation",	"Extra")
#############################################################
### VEP
Annotated <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/3-Annotation/2-Fichiers_annotEs/echantillon_fusionnE-annotE.txt",
                        header = FALSE, sep = "\t")
colnames( Annotated ) <- c( "Uploaded_variation",	"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence",	"cDNA_position",
                            "CDS_position",	"Protein_position",	"Amino_acids",	"Codons",	"Existing_variation",	"Extra" )
Annotated["sampleID"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[5] } ) )
Annotated["Patient"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[6] } ) )
Annotated["Batch"] <- ifelse( (substr(Annotated$sampleID, 1, 3) == "D21"), "batch2", "batch1" )
Annotated["Chr"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[1] } ) )
Annotated["pos"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[2] } ) )
Annotated["ref"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[3] } ) )
Annotated["alt"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[4] } ) )

library(comprehenr)
#synonymous_variant <- c( "synonymous_variant", "start_retained_variant", "stop_retained_variant" )
synonymous_variant <- c("synonymous_variant", "non_coding_transcript_variant", "start_retained_variant", "non_coding_transcript_exon_variant", "stop_retained_variant")
#nonsynonymous_variant <- c("missense_variant", "redundant_inserted_stop_gained", "start_lost", "stop_gained", "stop_lost")
#nonsynonymous_variant <- c( "nonsynonymous_variant", "missense_variant", "conservative_missense_variant", "non_conservative_missense_variant", "rare_amino_acid_variant", "pyrrolysine_loss", "selenocysteine_loss", "redundant_inserted_stop_gained", "start_lost", "stop_gained", "stop_gained_NMD_escaping", "stop_gained_NMD_triggering", "stop_lost" )
nonsynonymous_variant <- c("nonsynonymous_variant", "conservative_missense_variant", "non_conservative_missense_variant", "stop_gained_NMD_escaping",
                           "stop_gained_NMD_triggering", "rare_amino_acid_variant", "pyrrolysine_loss", "selenocysteine_loss", "redundant_inserted_stop_gained",
                           "frameshift_variant", "inframe_deletion", "splice_acceptor_variant", "stop_lost", "incomplete_terminal_codon_variant", "inframe_insertion",
                           "missense_variant", "protein_altering_variant", "splice_donor_variant", "start_lost", "stop_gained")
Annotated["Synonyme"] <- unlist( lapply( 1:nrow(Annotated), function(i) { ifelse( sum( to_vec(for(j in 1:length( unlist( strsplit( Annotated$Consequence[i], "," ) ) )) unlist( strsplit( Annotated$Consequence[i], "," ) )[j] %in% synonymous_variant ) ) >=1, 1, 0 ) } ) )
Annotated["Nonsynonyme"] <- unlist( lapply( 1:nrow(Annotated), function(i) { ifelse( sum( to_vec(for(j in 1:length( unlist( strsplit( Annotated$Consequence[i], "," ) ) )) unlist( strsplit( Annotated$Consequence[i], "," ) )[j] %in% nonsynonymous_variant ) ) >=1, 1, 0 ) } ) )

Annotated["impact"] <- unlist(lapply( 1:nrow(Annotated), function(i) { ifelse( Annotated$Nonsynonyme[i]==1, "Nonsynonymous", ifelse( Annotated$Synonyme[i]==1, "Synonymous", "Other" ) ) } ))


### dndscv
# Annotated <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/6-Ratio_dNdS/1-Sortie/mutations_annotation_by_dndscv.txt",
#                                header = TRUE, sep = "\t")
# Annotated$gene <- unlist( sapply( 1:nrow(Annotated), function(i) { unlist(strsplit(Annotated$gene[i], ":"))[2] } ) )
# Annotated["Patient"] <- unlist( lapply( 1:nrow(Annotated), function(i) { patients$numero[grep(Annotated$sampleID[i], patients$IDs)] } ) )
###

sample_ <- data.frame(id=NA, patient=NA, age=NA)
sample_ <- sample_[0,]
for (i in 1:nrow(patients)) {
  lines_ <- unlist(strsplit(patients$IDs[i], ","))
  for (j in 1:length(lines_)) {
    sample_[nrow(sample_)+1,] <- c(lines_[j], patients$numero[i], patients$Age[i])
  }
}

sample_["all_mutations"] <- unlist( lapply( 1:nrow(sample_), function(i) { nrow( subset( Annotated, sampleID==sample_$id[i] ) ) } ) )
sample_["Nonsynonymous"] <- unlist( lapply( 1:nrow(sample_), function(i) { nrow( subset( Annotated, (sampleID==sample_$id[i] & impact=="Nonsynonymous") ) ) } ) )
sample_["Synonymous"] <- unlist( lapply( 1:nrow(sample_), function(i) { nrow( subset( Annotated, (sampleID==sample_$id[i] & impact=="Synonymous") ) ) } ) )

sample_$age <- as.numeric(sample_$age)
sample_$age[15] <- sample_$age[15]+0.5
sample_$age[17] <- sample_$age[17]+0.5
sample_$age[18] <- sample_$age[18]+1.25
sample_$age[21] <- sample_$age[21]+0.5
sample_$age[22] <- sample_$age[22]+1

#write.table( patients, file="./1-Sortie/patients_updated.txt", col.names=T, quote=F, row.names=F, sep = "\t")

############################################################
sample_ <- transform( sample_, Nonsyn_rate = Nonsynonymous/all_mutations )

cov <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/0-Metadata/cov.txt", header = TRUE,
                  sep = "\t")
sample_["Cov_mean"] <- unlist( lapply( 1:nrow(sample_), function(i) { cov$Cov_mean[which(cov$Sample==sample_$id[i])] } ) )
sample_["N_norm"] <- sample_$all_mutations/sample_$Cov_mean
sample_["Nonsyn_norm"] <- sample_$Nonsynonymous/sample_$Cov_mean

sample_ordered <- sample_[order(sample_$age),]
sample_ordered$age <- as.numeric(sample_ordered$age)
sample_ordered$all_mutations <- as.numeric(sample_ordered$all_mutations)
sample_ordered$Nonsynonymous <- as.numeric(sample_ordered$Nonsynonymous)
sample_ordered$Synonymous <- as.numeric(sample_ordered$Synonymous)


#############################################################
three_patients <- sample_ordered[sample_ordered$patient %in% c(7, 8, 10),]

mutations <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/2-Filtre/2-Fichiers_filtrEs_avec_toutes_les_colonnes/echantillon_fusionnE-filtrE.txt",
                        header = FALSE, sep = "\t")
colnames(mutations) <- c( "CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "Details", "AF", "Sample", "Patient" )

mutations["id"] <- paste(mutations$CHROM, mutations$POS, mutations$REF, mutations$ALT, sep = "_")
three_patients <- mutations[mutations$Patient %in% c(7, 8, 10),]

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

## Entre échantillon du patient 8
patient_8 <- mutations[mutations$Patient == 8,]
D210288 <- patient_8[patient_8$Sample == "D210288",]$id
D210289 <- patient_8[patient_8$Sample == "D210289",]$id
D210290 <- patient_8[patient_8$Sample == "D210290",]$id

jaccard(D210288, D210289)
jaccard(D210289, D210290)

x <- list(
  Echantillon1 = D210288, 
  Echantillon2 = D210289
)
y <- list(
  Echantillon2 = D210289,
  Echantillon3 = D210290
)
library(ggvenn)
png(file =  "./2-Images/1-Patient8_sample1-2.png", width = 900, height = 800)
ggvenn(
  x,
  text_size = 10,
  auto_scale = TRUE,
  digits = 1,
  fill_color = c("red", "blue"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title = "Prélèvements espacés de 6 mois (Cytobrosse)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
dev.off()

png(file =  "./2-Images/2-Patient8_sample2-3.png", width = 900, height = 800)
ggvenn(
  y, 
  text_size = 10,
  auto_scale = TRUE,
  digits = 1,
  fill_color = c("blue", "green"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title = "Prélèvements espacés de 9 mois (Cytobrosse)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
dev.off()

## Entre échantillon du patient 7
patient_7 <- mutations[mutations$Patient == 7,]
D210285 <- patient_7[patient_7$Sample == "D210285",]$id
D210327 <- patient_7[patient_7$Sample == "D210327",]$id
D210328 <- patient_7[patient_7$Sample == "D210328",]$id

x <- list(
  Echantillon1 = D210285, 
  Echantillon2 = D210327
)
y <- list(
  Echantillon2 = D210327,
  Echantillon3 = D210328
)
library(ggvenn)
png(file =  "./2-Images/3-Patient7_sample1-2.png", width = 900, height = 800)
ggvenn(
  x, 
  text_size = 10,
  auto_scale = TRUE,
  digits = 1,
  fill_color = c("red", "blue"),
) +
  labs(title = "Prélèvements au même moment (Cytobrosse)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
dev.off()
png(file =  "./2-Images/4-Patient7_sample2-3.png", width = 900, height = 800)
ggvenn(
  y, 
  text_size = 10,
  auto_scale = TRUE,
  digits = 1,
  fill_color = c("blue", "green"),
) +
  labs(title = "Prélèvements espacés de 6 mois (Cytobrosse)") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
dev.off()

############################################################
## Brosse-brosse interpatient
Jaccard_brosse_brosse_interpatient <- data.frame(id1=character(), Jaccard_index=numeric(), id2=character())
p <- c("D181208", "D181210", "D181211", "D181212", "D181213", "D181215", "D210285",
       "D210288", "D210297", "D210294", "D210332", "D210334", "D210335", "D210338", "D210339")
k = 1
done <- c()
for (i in 1:15) {
  for (j in 1:15) {
    if (i!=j) {
      m <- c(p[i], p[j])
      m_sorted <- sort(m, decreasing = FALSE)
      m_pasted <- paste(m_sorted[1], m_sorted[2], sep = "")
      if (! m_pasted %in% done) {
        p1 <- subset(mutations, Sample==p[i])
        p2 <- subset(mutations, Sample==p[j])
        Jaccard_brosse_brosse_interpatient[k,] <- c(p[i], jaccard(p1$id, p2$id), p[j])
        k = k+1
        done <- c(done, m_pasted)
      }
    }
  }
}
Jaccard_brosse_brosse_interpatient["VS"] <- "Brosse VS Brosse (patients différents)"
# Jaccard_brosse_brosse_interpatient["VS"] <- ifelse(substr(Jaccard_brosse_brosse_interpatient$id1, 1, 3)=="D18", "batch1_VS_", "batch2_VS")
# Jaccard_brosse_brosse_interpatient$VS <- ifelse(substr(Jaccard_brosse_brosse_interpatient$id2, 1, 3)=="D18",
#                                                paste(Jaccard_brosse_brosse_interpatient$VS, "batch1", sep=""), 
#                                 paste(Jaccard_brosse_brosse_interpatient$VS, "batch2", sep=""))
# Jaccard_brosse_brosse_interpatient$Jaccard_index <- as.numeric(Jaccard_brosse_brosse_interpatient$Jaccard_index)

## Brosse-brosse même patient
Jaccard_brosse_brosse_memepatient <- data.frame(id1=character(), Jaccard_index=numeric(), id2=character())
p <- list(c("D210285","D210327","D210328"), c("D210288","D210289","D210290"), c("D210294","D210295","D210296"))
k = 1
for (i in 1:length(p)) {
  p1 <- subset(mutations, Sample==p[[i]][1])
  p2 <- subset(mutations, Sample==p[[i]][2])
  p3 <- subset(mutations, Sample==p[[i]][3])
  Jaccard_brosse_brosse_memepatient[k,] <- c(p[[i]][1], jaccard(p1$id, p2$id), p[[i]][2])
  k = k+1
  Jaccard_brosse_brosse_memepatient[k,] <- c(p[[i]][1], jaccard(p1$id, p3$id), p[[i]][3])
  k = k+1
  Jaccard_brosse_brosse_memepatient[k,] <- c(p[[i]][2], jaccard(p2$id, p3$id), p[[i]][3])
  k = k+1
}
Jaccard_brosse_brosse_memepatient["VS"] <- "Brosse VS Brosse (même patient)"

## Biopsie-brosse même patient
Jaccard_biopsie_brosse_memepatient <- data.frame(id1=character(), Jaccard_index=numeric(), id2=character())
p <- list(c("D181208","D202062"), c("D181210","D202061"), c("D181211","D202063"),
          c("D181212","D202064"), c("D181213","D202065"), c("D181215","D202066"))
k = 1
for (i in 1:length(p)) {
  p1 <- subset(mutations, Sample==p[[i]][1])
  p2 <- subset(mutations, Sample==p[[i]][2])
  Jaccard_biopsie_brosse_memepatient[k,] <- c(p[[i]][1], jaccard(p1$id, p2$id), p[[i]][2])
  k = k+1
}
Jaccard_biopsie_brosse_memepatient["VS"] <- "Biopsie VS Brosse (même patient)"

## TESTS
test_1 <- wilcox.test(as.numeric(Jaccard_biopsie_brosse_memepatient$Jaccard_index),
                      as.numeric(Jaccard_brosse_brosse_memepatient$Jaccard_index), paired = FALSE)
p1_ <- test_1$p.value

test_2 <- wilcox.test(as.numeric(Jaccard_biopsie_brosse_memepatient$Jaccard_index),
                      as.numeric(Jaccard_brosse_brosse_interpatient$Jaccard_index), paired = FALSE)
p2_ <- test_2$p.value

test_3 <- wilcox.test(as.numeric(Jaccard_brosse_brosse_memepatient$Jaccard_index),
                      as.numeric(Jaccard_brosse_brosse_interpatient$Jaccard_index), paired = FALSE)
p3 <- test_3$p.value

##
Jaccard_table <- rbind(Jaccard_brosse_brosse_interpatient, Jaccard_brosse_brosse_memepatient, Jaccard_biopsie_brosse_memepatient)
Jaccard_table$Jaccard_index <- as.numeric(Jaccard_table$Jaccard_index)

## plot
library(tidyverse)
library(hrbrthemes)
library(viridis)
png(file =  "./2-Images/6-Boxplot_indice_Jaccard.png", width = 900, height = 800)
Jaccard_table %>%
  ggplot( aes(x=VS, y=Jaccard_index, color = VS)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), cex=3) +
  # geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0, seed=1234),
  #            pch=21,
  #            aes(fill=factor(VS), size = 4),
  #            show.legend = F) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    legend.title=element_blank(),
    legend.text=element_text(size=20),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 20),
  ) +
  ylab("Indice de Jaccard") +
  xlab("Echantillons croisés")
dev.off()


############################################################

library(stats)
sample_ordered["Batch"] <- ifelse(substr(sample_ordered$id, 1, 3)=="D21", "batch1",
                                  "batch2")
reg <- lm(sample_ordered$N_norm~sample_ordered$age)
sample_ordered$Batch <- as.factor(sample_ordered$Batch)
png(file =  "./2-Images/5-Mutation_selon_age.png", width = 900, height = 800)
par(mar=c(5, 6, 5, 1))
plot(sample_ordered$N_norm~sample_ordered$age, col=sample_ordered$Batch, pch=19,
     xlab="Age", ylab="Nombre de mutation total normalisé \n par la couverture moyenne pour chaque échantillon",
     cex.lab=1.5)
abline(reg, col="blue")
legend("topright",legend=levels(sample_ordered$Batch),col=1:2,pch=19, cex=2)
dev.off()

summary(reg)
cor(sample_ordered$N_norm, sample_ordered$age)
