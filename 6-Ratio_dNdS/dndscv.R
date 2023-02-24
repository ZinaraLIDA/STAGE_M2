rm(list = ls())

########################################################################################################################################################################
setwd("~/ZINARA/Gnomic/1-StageM2/6-Ratio_dNdS")

# Importation des jeux de données
dataframe_all_sample_filtered <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/5-Stats/1-text/Apres_filtre_avec_toutes_les_colonnes.txt",
                                                 header = TRUE, stringsAsFactors = FALSE, sep = "\t"  )
# dataframe_all_sample_filtered <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/5-Stats/1-text/Avant_filtre_avec_toutes_les_colonnes.txt", # Non filtré
#                                              header = TRUE, stringsAsFactors = FALSE, sep = "\t"  )
# dataframe_all_sample_filtered <- subset(dataframe_all_sample_filtered, Batch == "batch2")
###################################
library("dndscv")
data_ <- dataframe_all_sample_filtered[, c( "Sample", "CHROM", "POS", "REF", "ALT" )]
colnames(data_) <- c( "sampleID", "chr", "pos", "ref", "mut" )

cds <- read.table( "./mart_export.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE )
cds[, "Chromosome.scaffold.name"] <- paste( "chr", unlist(cds["Chromosome.scaffold.name"]), sep = "")
write.table(cds, file="cds.txt", col.names=T, quote=F, row.names=F, sep = "\t")

#buildref(cdsfile="./cds.txt", genomefile="./hg38-au.pri.fa", outfile = "refcds.rda", excludechrs="MT", useids = TRUE)

#data_$chr = gsub("chr","",as.vector(data_$chr))
# library(dplyr)
# data_ <- distinct(data_, chr, pos, ref, mut, .keep_all= TRUE)
dndsout = dndscv(data_, refdb = "./refcds.rda",
                 max_coding_muts_per_sample = Inf,
                 max_muts_per_gene_per_sample = Inf,
                 cv=NULL
                 #,gene_list = as.character(gene_list$V1)
                 )
sel_cv = dndsout$sel_cv
sel_cv["gene_name"] <- unlist( sapply( 1:nrow(sel_cv), function(i) { unlist(strsplit(sel_cv$gene_name[i], ":"))[2] } ) )

signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

#########################################################################
drivers_genes_list <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/0-Metadata/IntOGen-DriverGenes_HNSC_2022.tsv", header = TRUE, sep = "\t", dec = ".",
                                  stringsAsFactors = FALSE )
most500_mutated_gene_list <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/5-Stats/1-text/500_genes_les_plus_mutes_parmi_les_genes_exprimes.txt",
                                         header = TRUE, sep = "\t", stringsAsFactors = FALSE )
most1000_mutated_unexpressed_gene_list <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/5-Stats/1-text/1000_genes_les_plus_mutes_parmi_les_genes_non_exprimes.txt",
                                           header = TRUE, sep = "\t", stringsAsFactors = FALSE )
###################################
data_drivers_genes <- sel_cv[sel_cv$gene_name %in% drivers_genes_list$Symbol,]
data_most500_mutated_genes <- sel_cv[sel_cv$gene_name %in% most500_mutated_gene_list$Genes,]
data_most1000_mutated_unexpressed_genes <- sel_cv[sel_cv$gene_name %in% most1000_mutated_unexpressed_gene_list$V1,]
###################################
write.table(sel_cv, file="./1-Sortie/sel_cv.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table(data_drivers_genes, file="./1-Sortie/data_drivers_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table(data_most500_mutated_genes, file="./1-Sortie/data_most500_mutated_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table(data_most1000_mutated_unexpressed_genes, file="./1-Sortie/data_most1000_mutated_unexpressed_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")
write.table(dndsout[["annotmuts"]], file="./1-Sortie/mutations_annotation_by_dndscv.txt", col.names=T,
            quote=F, row.names=F, sep = "\t")
write.table(known_cancergenes, file="./1-Sortie/known_cancergenes.txt", col.names=F,
            quote=F, row.names=F, sep = "\t")
#####################################
# Les gènes avec un avantage sélectif positif
write.table(signif_genes, file="./2-Significatif/signif_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

signif_drivers_genes <- signif_genes[signif_genes$gene_name %in% drivers_genes_list$Symbol,]
write.table(signif_drivers_genes, file="./2-Significatif/signif_drivers_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

signif_most500_mutated_genes <- signif_genes[signif_genes$gene_name %in% most500_mutated_gene_list$Genes,]
write.table(signif_most500_mutated_genes, file="./2-Significatif/signif_most500_mutated_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

signif_most1000_mutated_unexpressed_genes <- signif_genes[signif_genes$gene_name %in% most1000_mutated_unexpressed_gene_list$Genes,]
write.table(signif_most1000_mutated_unexpressed_genes, file="./2-Significatif/signif_most1000_mutated_unexpressed_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")
#####################################
# Les gènes sous sélection négative
dndsout$sel_cv$pmis_negsel = dndsout$sel_cv$pmis_cv / 2 # Initialising the one-sided vector of p-values for missense mutations
dndsout$sel_cv$pmis_negsel[dndsout$sel_cv$wmis_cv>1] = 1 # Setting one sided p-values for genes with dN/dS>1 to 1
dndsout$sel_cv$qmis_negsel = p.adjust(dndsout$sel_cv$pmis_negsel, method="BH") # FDR adjustment (q-values)
head(dndsout$sel_cv[order(dndsout$sel_cv$qmis_negsel),]) # Looking at the most significant genes for negative selection

most1000_negative_selected <- dndsout$sel_cv[order(dndsout$sel_cv$qmis_negsel),][1:1000,]
most1000_negative_selected$gene_name<- unlist( sapply( 1:nrow(most1000_negative_selected),
                                                       function(i) { unlist(strsplit(most1000_negative_selected$gene_name[i], ":"))[2] } ) )

write.table(most1000_negative_selected, file="./3-Selection_negative/most1000_negative_selected_genes.txt", col.names=T,
            quote=F, row.names=F, sep = "\t")
######################################


x <- list(
  'Drivers CECO' = drivers_genes_list$Symbol, 
  'Gènes significatifs' = signif_genes$gene_name
)

y <- list(
  'Gènes du cancer connus' = known_cancergenes, 
  'Gènes significatifs' = signif_genes$gene_name
)

z <- list(
  "Gènes les plus mutés parmi ceux exprimés dans l'épithélium de la muqueuse orale" = data_most500_mutated_genes$gene_name, 
  'Gènes significatifs' = signif_genes$gene_name
)

library(ggvenn)
png(file =  "./4-Images/1-Entre_drivers_CECO_et_genes_signif.png", width = 900, height = 800)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title = "Gènes en commun entre les gènes drivers CECO et les gènes significatifs") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

png(file =  "./4-Images/2-Entre_genes_du_cancer_et_genes_signif.png", width = 900, height = 800)
ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title = "Gènes en commun entre les gènes connus du cancer et les gènes significatifs") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

png(file =  "./4-Images/3-Entre_muqueuse_orale_et_genes_signif.png", width = 900, height = 800)
ggvenn(
  z, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title = "Gènes en commun entre les gènes exprimés \n dans l'épithélium de la muqueuse orale et les gènes significatifs") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()


######################################





data_new <- data_
data_new["id"] <- paste(data_$sampleID, data_$chr, data_$pos, data_$ref, data_$mut, sep = "_")

annot_mut <- dndsout[["annotmuts"]]
annot_mut["id"] <- paste(annot_mut$sampleID, annot_mut$chr, annot_mut$pos, annot_mut$ref, annot_mut$mut, sep = "_")
annot_mut$pos <- as.numeric(annot_mut$pos)

diff <- data_new[!(data_new$id %in% annot_mut$id),]

all_ <- dataframe_all_sample_filtered
all_["id"] <- paste(all_$Sample, all_$CHROM, all_$POS, all_$REF, all_$ALT, sep = "_")
diff_all <- all_[(all_$id %in% diff$id),]





##
table(diff_all$CHROM)
barplot(table(diff_all$CHROM))
hist(diff_all$QUAL)
quantile(diff_all$QUAL) 
hist(diff_all$AF)
barplot(table(diff_all$Sample))
hist(diff_all$DP)
table(diff_all$Batch)
quantile(diff_all$AF)
quantile(diff_all$AO)
mean(diff_all$QUAL)
mean(diff_all$AF)