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
# Annotated <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/3-Annotation/2-Fichiers_annotEs/echantillon_fusionnE-annotE.txt",
#                         header = FALSE, sep = "\t")
# colnames( Annotated ) <- c( "Uploaded_variation",	"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence",	"cDNA_position",
#                             "CDS_position",	"Protein_position",	"Amino_acids",	"Codons",	"Existing_variation",	"Extra" )
# Annotated["Sample"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[5] } ) )
# Annotated["Patient"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[6] } ) )
# Annotated["Batch"] <- ifelse( (substr(Annotated$Sample, 1, 3) == "D21"), "batch2", "batch1" )
# Annotated["Chr"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[1] } ) )
# Annotated["pos"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[2] } ) )
# Annotated["ref"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[3] } ) )
# Annotated["alt"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[4] } ) )
# library(dplyr)
# Annotated <- distinct(Annotated, Chr, pos, ref, alt, Feature, .keep_all= TRUE)
# 
# library(comprehenr)
# #synonymous_variant <- c( "synonymous_variant", "start_retained_variant", "stop_retained_variant" )
# synonymous_variant <- c("synonymous_variant", "non_coding_transcript_variant", "start_retained_variant", "non_coding_transcript_exon_variant", "stop_retained_variant")
# #nonsynonymous_variant <- c("missense_variant", "redundant_inserted_stop_gained", "start_lost", "stop_gained", "stop_lost")
# #nonsynonymous_variant <- c( "nonsynonymous_variant", "missense_variant", "conservative_missense_variant", "non_conservative_missense_variant", "rare_amino_acid_variant", "pyrrolysine_loss", "selenocysteine_loss", "redundant_inserted_stop_gained", "start_lost", "stop_gained", "stop_gained_NMD_escaping", "stop_gained_NMD_triggering", "stop_lost" )
# nonsynonymous_variant <- c("nonsynonymous_variant", "conservative_missense_variant", "non_conservative_missense_variant", "stop_gained_NMD_escaping", "stop_gained_NMD_triggering", "rare_amino_acid_variant", "pyrrolysine_loss", "selenocysteine_loss", "redundant_inserted_stop_gained", "frameshift_variant", "inframe_deletion", "splice_acceptor_variant", "stop_lost", "incomplete_terminal_codon_variant", "inframe_insertion", "missense_variant", "protein_altering_variant", "splice_donor_variant", "start_lost", "stop_gained")
# Annotated["Synonyme"] <- unlist( lapply( 1:nrow(Annotated), function(i) { ifelse( sum( to_vec(for(j in 1:length( unlist( strsplit( Annotated$Consequence[i], "," ) ) )) unlist( strsplit( Annotated$Consequence[i], "," ) )[j] %in% synonymous_variant ) ) >=1, 1, 0 ) } ) )
# Annotated["Nonsynonyme"] <- unlist( lapply( 1:nrow(Annotated), function(i) { ifelse( sum( to_vec(for(j in 1:length( unlist( strsplit( Annotated$Consequence[i], "," ) ) )) unlist( strsplit( Annotated$Consequence[i], "," ) )[j] %in% nonsynonymous_variant ) ) >=1, 1, 0 ) } ) )
# 
# Annotated["impact"] <- lapply( 1:nrow(Annotated), function(i) { ifelse( Annotated$Nonsynonyme[i]==1, "Nonsynonymous", ifelse( Annotated$Synonyme[i]==1, "Synonymous", "Other" ) ) } )


### dndscv
Annotated <- read.table(file = "~/ZINARA/Gnomic/1-StageM2/6-Ratio_dNdS/1-Sortie/mutations_annotation_by_dndscv.txt",
                               header = TRUE, sep = "\t")
Annotated$gene <- unlist( sapply( 1:nrow(Annotated), function(i) { unlist(strsplit(Annotated$gene[i], ":"))[2] } ) )

library(dplyr)
Annotated <- distinct(Annotated, chr, pos, ref, mut, .keep_all= TRUE)

Annotated["Patient"] <- unlist( lapply( 1:nrow(Annotated), function(i) { patients$numero[grep(Annotated$sampleID[i], patients$IDs)] } ) )
###


patients["all_mutations"] <- unlist( lapply( 1:nrow(patients), function(i) { nrow( subset( Annotated, Patient==patients$numero[i] ) ) } ) )
patients["Nonsynonymous"] <- unlist( lapply( 1:nrow(patients), function(i) { nrow( subset( Annotated, (Patient==patients$numero[i] & impact!="Synonymous") ) ) } ) )
patients["Synonymous"] <- unlist( lapply( 1:nrow(patients), function(i) { nrow( subset( Annotated, (Patient==patients$numero[i] & impact=="Synonymous") ) ) } ) )

write.table( patients, file="./1-Sortie/patients_updated.txt", col.names=T, quote=F, row.names=F, sep = "\t")

############################################################
patients <- transform( patients, Nonsyn_rate = Nonsynonymous/all_mutations ) 
patients_ordered <- patients[order(patients$Age),]
patients_ordered$Age <- as.numeric(patients_ordered$Age)
patients_ordered$all_mutations <- as.numeric(patients_ordered$all_mutations)
patients_ordered$Nonsynonymous <- as.numeric(patients_ordered$Nonsynonymous)
patients_ordered$Synonymous <- as.numeric(patients_ordered$Synonymous)

library(stats)
reg <- lm(patients_ordered$Nonsyn_rate~patients_ordered$Age)
plot(patients_ordered$Nonsyn_rate~patients_ordered$Age)
abline(reg, col="red")
summary(reg)
