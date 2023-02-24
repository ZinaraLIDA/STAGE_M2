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
Annotated["Sample"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[5] } ) )
Annotated["Patient"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[6] } ) )
Annotated["Batch"] <- ifelse( (substr(Annotated$Sample, 1, 3) == "D21"), "batch2", "batch1" )
Annotated["Chr"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[1] } ) )
Annotated["pos"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[2] } ) )
Annotated["ref"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[3] } ) )
Annotated["alt"] <- unlist( lapply( 1:nrow(Annotated), function(i) { unlist( strsplit( Annotated$Uploaded_variation[i], "_" ) )[4] } ) )
# library(dplyr)
# Annotated <- distinct(Annotated, Chr, pos, ref, alt, Feature, .keep_all= TRUE)

library(comprehenr)
#---synonymous_variant <- c( "synonymous_variant", "start_retained_variant", "stop_retained_variant" )
synonymous_variant <- c("synonymous_variant", "non_coding_transcript_variant", "start_retained_variant", "non_coding_transcript_exon_variant", "stop_retained_variant")
#---nonsynonymous_variant <- c("missense_variant", "redundant_inserted_stop_gained", "start_lost", "stop_gained", "stop_lost")
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
# 
# library(dplyr)
# Annotated <- distinct(Annotated, chr, pos, ref, mut, .keep_all= TRUE)
# 
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

sample_["all_mutations"] <- unlist( lapply( 1:nrow(sample_), function(i) { nrow( subset( Annotated, Sample==sample_$id[i] ) ) } ) )
sample_["Nonsynonymous"] <- unlist( lapply( 1:nrow(sample_), function(i) { nrow( subset( Annotated, (Sample==sample_$id[i] & impact=="Nonsynonymous") ) ) } ) )
sample_["Synonymous"] <- unlist( lapply( 1:nrow(sample_), function(i) { nrow( subset( Annotated, (Sample==sample_$id[i] & impact=="Synonymous") ) ) } ) )

sample_$age <- as.numeric(sample_$age)
sample_$age[15] <- sample_$age[15]+0.5
sample_$age[17] <- sample_$age[17]+0.5
sample_$age[18] <- sample_$age[18]+1.25
sample_$age[21] <- sample_$age[21]+0.5
sample_$age[22] <- sample_$age[22]+1

#write.table( patients, file="./1-Sortie/patients_updated.txt", col.names=T, quote=F, row.names=F, sep = "\t")

############################################################
sample_ <- transform( sample_, Nonsyn_rate = Nonsynonymous/all_mutations ) 
sample_ordered <- sample_[order(sample_$age),]
sample_ordered$age <- as.numeric(sample_ordered$age)
sample_ordered$all_mutations <- as.numeric(sample_ordered$all_mutations)
sample_ordered$Nonsynonymous <- as.numeric(sample_ordered$Nonsynonymous)
sample_ordered$Synonymous <- as.numeric(sample_ordered$Synonymous)

library(stats)
my_col <- c("red","blue")
sample_ordered["Batch"] <- ifelse(substr(sample_ordered$id, 1, 3)=="D21", "batch1",
                                  "batch2")
reg <- lm(sample_ordered$all_mutations~sample_ordered$age)
plot(sample_ordered$all_mutations~sample_ordered$age,
     col = my_col[sample_ordered$Batch])
abline(reg, col="red")
summary(reg)
