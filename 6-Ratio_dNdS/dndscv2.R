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
write.table(signif_genes, file="./2-Significatif/signif_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

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
data_most1000_mutated_unexpressed_genes <- sel_cv[sel_cv$gene_name %in% most1000_mutated_unexpressed_gene_list$Genes,]
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
signif_drivers_genes <- signif_genes[signif_genes$gene_name %in% drivers_genes_list$Symbol,]
write.table(signif_drivers_genes, file="./2-Significatif/signif_drivers_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

signif_most500_mutated_genes <- signif_genes[signif_genes$gene_name %in% most500_mutated_gene_list$Genes,]
write.table(signif_most500_mutated_genes, file="./2-Significatif/Significatifs_et_exprimes.txt", col.names=T, quote=F, row.names=F, sep = "\t")

signif_most1000_mutated_unexpressed_genes <- signif_genes[signif_genes$gene_name %in% most1000_mutated_unexpressed_gene_list$Genes,]
write.table(signif_most1000_mutated_unexpressed_genes, file="./2-Significatif/Significatifs_et_nonexprimes.txt", col.names=T, quote=F, row.names=F, sep = "\t")
#####################################
# Les gènes avec un avantage sélectif positif
pos_signif_genes = sel_cv[sel_cv$qglobal_cv<0.1 & sel_cv$wmis_cv>1, c("gene_name","qglobal_cv")]
write.table(pos_signif_genes, file="./2-Significatif/pos_signif_genes.txt", col.names=T, quote=F, row.names=F, sep = "\t")
pos_and_driversGenes <- drivers_genes_list[drivers_genes_list$Symbol %in% pos_signif_genes$gene_name,]$Symbol
write.table(pos_and_driversGenes, file="./2-Significatif/Selection_positive_et_Drivers.txt", col.names=F, quote=F, row.names=F, sep = "\t")
pos_and_expressedGenes <- pos_signif_genes[pos_signif_genes$gene_name %in% most500_mutated_gene_list$Genes,]$gene_name
pos_and_nonexpressedGenes <- pos_signif_genes[pos_signif_genes$gene_name %in% most1000_mutated_unexpressed_gene_list$Genes,]$gene_name
write.table(pos_and_expressedGenes, file="./2-Significatif/Selection_positive_et_exprimes.txt", col.names=F, quote=F, row.names=F, sep = "\t")
write.table(pos_and_nonexpressedGenes, file="./2-Significatif/Selection_positive_et_nonexprimes.txt", col.names=F, quote=F, row.names=F, sep = "\t")
#####################################
# Les gènes sous sélection négative
neg_signif_genes = sel_cv[sel_cv$qglobal_cv<0.1 & sel_cv$wmis_cv<1, c("gene_name","qglobal_cv")]
write.table(neg_signif_genes, file="./2-Significatif/Selection_negative.txt", col.names=T, quote=F, row.names=F, sep = "\t")
neg_and_driversGenes <- drivers_genes_list[drivers_genes_list$Symbol %in% neg_signif_genes$gene_name,]$Symbol
write.table(neg_and_driversGenes, file="./2-Significatif/Selection_negative_et_Drivers.txt", col.names=F, quote=F, row.names=F, sep = "\t")
neg_and_expressedGenes <- neg_signif_genes[neg_signif_genes$gene_name %in% most500_mutated_gene_list$Genes,]$gene_name
neg_and_nonexpressedGenes <- neg_signif_genes[neg_signif_genes$gene_name %in% most1000_mutated_unexpressed_gene_list$Genes,]$gene_name
write.table(neg_and_expressedGenes, file="./2-Significatif/Selection_negative_et_exprimes.txt", col.names=F, quote=F, row.names=F, sep = "\t")
write.table(neg_and_nonexpressedGenes, file="./2-Significatif/Selection_negative_et_nonexprimes.txt", col.names=F, quote=F, row.names=F, sep = "\t")
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
png(file =  "./3-Images/1-Entre_drivers_CECO_et_genes_signif.png", width = 900, height = 800)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title = "Gènes en commun entre les gènes drivers CECO et les gènes significatifs") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

png(file =  "./3-Images/2-Entre_genes_du_cancer_et_genes_signif.png", width = 900, height = 800)
ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title = "Gènes en commun entre les gènes connus du cancer et les gènes significatifs") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()

png(file =  "./3-Images/3-Entre_muqueuse_orale_et_genes_signif.png", width = 900, height = 800)
ggvenn(
  z, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title = "Gènes en commun entre les gènes exprimés \n dans l'épithélium de la muqueuse orale et les gènes significatifs") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dev.off()


#########################################################################################################################################################################
#########################################################################################################################################################################
############################################################################ AVEC CORRECTION ############################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################

# 62 GENES DRIVERS
##################

data_drivers_genes_corrected <- data_drivers_genes
library(multtest)
matrix_ <- mt.rawp2adjp(data_drivers_genes_corrected$pglobal_cv, proc = "BH", alpha = 0.05)
adj <- matrix_$adjp[order(matrix_$index),][,2]
data_drivers_genes_corrected["pglobalcv_adjusted"] <- adj

data_drivers_genes_corrected <- transform(data_drivers_genes_corrected, n_mut = n_syn+n_mis+n_non+n_spl+n_ind)
data_drivers_genes_corrected <- data_drivers_genes_corrected[, c("gene_name", "n_mut", "wmis_cv", "wnon_cv", "pglobalcv_adjusted")]

data_drivers_genes_corrected2 <- subset(data_drivers_genes_corrected, n_mut>=5 & pglobalcv_adjusted<0.05)
colnames(data_drivers_genes_corrected2) <- c("Gène", "nombre_de_mutations", "dNdS_faux_sens", "dNdS_non_sens", "p_globale_ajustée")

library(formattable)
widget_formattable = formattable(data_drivers_genes_corrected2, list(
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens>1)) ~ color_tile("salmon","red")
))
widget_formattable

html_header="
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
</head>
<body>
"
html_table = format_table(data_drivers_genes_corrected2, list(
  area(col = 3, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red")
))
write(paste(html_header, html_table, sep=""), "./4-Correction/62/Supplementary_table_2_CrypticMAPSs.html")
library(htmltools)
library(webshot)
webshot::install_phantomjs()
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(widget_formattable, "./4-Correction/62/Supplementary_table_2_CrypticMAPSs.png")
write.table(data_drivers_genes_corrected2, file="./4-Correction/62/Genes_drivers2.csv", col.names=T, quote=F, row.names=F, sep = "\t")

#####

data_drivers_genes_corrected3 <- subset(data_drivers_genes_corrected, n_mut>=5 & pglobalcv_adjusted>=0.05)
colnames(data_drivers_genes_corrected3) <- c("Gène", "nombre_de_mutation", "dNdS_faux_sens", "dNdS_non_sens", "p_valeur_ajustée")

widget_formattable = formattable(data_drivers_genes_corrected3, list(
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens>1)) ~ color_tile("salmon","red")
))
widget_formattable

html_header="
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
</head>
<body>
"
html_table = format_table(data_drivers_genes_corrected3, list(
  area(col = 3, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red")
))
write(paste(html_header, html_table, sep=""), "./4-Correction/62/Supplementary_table_3_CrypticMAPSs.html")
library(htmltools)
library(webshot)
webshot::install_phantomjs()
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(widget_formattable, "./4-Correction/62/Supplementary_table_3_CrypticMAPSs.png")
write.table(data_drivers_genes_corrected3, file="./4-Correction/62/Genes_drivers3.csv", col.names=T, quote=F, row.names=F, sep = "\t")

# 500 GENES LES PLUS MUTES ET EXPRIMES
######################################

genes500_corrected <- data_most500_mutated_genes
library(multtest)
matrix_ <- mt.rawp2adjp(genes500_corrected$pglobal_cv, proc = "BH", alpha = 0.05)
adj <- matrix_$adjp[order(matrix_$index),][,2]
genes500_corrected["pglobalcv_adjusted"] <- adj

genes500_corrected <- transform(genes500_corrected, n_mut = n_syn+n_mis+n_non+n_spl+n_ind)
genes500_corrected <- genes500_corrected[, c("gene_name", "n_mut", "wmis_cv", "wnon_cv", "pglobalcv_adjusted")]

genes500_corrected2 <- subset(genes500_corrected, n_mut>=5 & pglobalcv_adjusted<0.05)
colnames(genes500_corrected2) <- c("Gène", "nombre_de_mutations", "dNdS_faux_sens", "dNdS_non_sens", "p_globale_ajustée")

library(formattable)
widget_formattable = formattable(genes500_corrected2, list(
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens>1)) ~ color_tile("salmon","red")
))
widget_formattable

html_header="
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
</head>
<body>
"
html_table = format_table(genes500_corrected2, list(
  area(col = 3, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red")
))
write(paste(html_header, html_table, sep=""), "./4-Correction/500/Supplementary_table_2_CrypticMAPSs.html")
library(htmltools)
library(webshot)
webshot::install_phantomjs()
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(widget_formattable, "./4-Correction/500/Supplementary_table_2_CrypticMAPSs.png")
write.table(genes500_corrected2, file="./4-Correction/500/Genes_drivers2.csv", col.names=T, quote=F, row.names=F, sep = "\t")

#####

genes500_corrected3 <- subset(genes500_corrected, n_mut>=5 & pglobalcv_adjusted>=0.05)
colnames(genes500_corrected3) <- c("Gène", "nombre_de_mutation", "dNdS_faux_sens", "dNdS_non_sens", "p_valeur_ajustée")

widget_formattable = formattable(genes500_corrected3, list(
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens>1)) ~ color_tile("salmon","red")
))
widget_formattable

html_header="
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
</head>
<body>
"
html_table = format_table(genes500_corrected3, list(
  area(col = 3, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red")
))
write(paste(html_header, html_table, sep=""), "./4-Correction/500/Supplementary_table_3_CrypticMAPSs.html")
library(htmltools)
library(webshot)
webshot::install_phantomjs()
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(widget_formattable, "./4-Correction/500/Supplementary_table_3_CrypticMAPSs.png")
write.table(genes500_corrected3, file="./4-Correction/500/Genes_drivers3.csv", col.names=T, quote=F, row.names=F, sep = "\t")
export_formattable(widget_formattable, "./4-Correction/500/Supplementary_table_3_CrypticMAPSs.xls")
write.table(widget_formattable, file="./4-Correction/500/Genes_drivers3.xls", col.names=T, quote=F, row.names=F, sep = "\t")

# 1000 GENES LES PLUS MUTES ET EXPRIMES
######################################

genes1000_corrected <- data_most1000_mutated_unexpressed_genes
library(multtest)
matrix_ <- mt.rawp2adjp(genes1000_corrected$pglobal_cv, proc = "BH", alpha = 0.05)
adj <- matrix_$adjp[order(matrix_$index),][,2]
genes1000_corrected["pglobalcv_adjusted"] <- adj

genes1000_corrected <- transform(genes1000_corrected, n_mut = n_syn+n_mis+n_non+n_spl+n_ind)
genes1000_corrected <- genes1000_corrected[, c("gene_name", "n_mut", "wmis_cv", "wnon_cv", "pglobalcv_adjusted")]

genes1000_corrected2 <- subset(genes1000_corrected, n_mut>=5 & pglobalcv_adjusted<0.05)
colnames(genes1000_corrected2) <- c("Gène", "nombre_de_mutations", "dNdS_faux_sens", "dNdS_non_sens", "p_globale_ajustée")

library(formattable)
widget_formattable = formattable(genes1000_corrected2, list(
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens>1)) ~ color_tile("salmon","red")
))
widget_formattable

html_header="
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
</head>
<body>
"
html_table = format_table(genes1000_corrected2, list(
  area(col = 3, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red")
))
write(paste(html_header, html_table, sep=""), "./4-Correction/1000/Supplementary_table_2_CrypticMAPSs.html")
library(htmltools)
library(webshot)
webshot::install_phantomjs()
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(widget_formattable, "./4-Correction/1000/Supplementary_table_2_CrypticMAPSs.png")
write.table(genes1000_corrected2, file="./4-Correction/1000/Genes_drivers2.csv", col.names=T, quote=F, row.names=F, sep = "\t")

#####

genes1000_corrected3 <- subset(genes1000_corrected, n_mut>=5 & pglobalcv_adjusted>=0.05)
colnames(genes1000_corrected3) <- c("Gène", "nombre_de_mutation", "dNdS_faux_sens", "dNdS_non_sens", "p_valeur_ajustée")

widget_formattable = formattable(genes1000_corrected3, list(
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$dNdS_faux_sens>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$dNdS_non_sens>1)) ~ color_tile("salmon","red")
))
widget_formattable

html_header="
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
</head>
<body>
"
html_table = format_table(genes1000_corrected3, list(
  area(col = 3, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 3, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red"),
  area(col = 4, row = which(widget_formattable$wmis_cv<=1)) ~ color_tile('blue', 'lightblue'),
  area(col = 4, row = which(widget_formattable$wmis_cv>1)) ~ color_tile("salmon","red")
))
write(paste(html_header, html_table, sep=""), "./4-Correction/1000/Supplementary_table_3_CrypticMAPSs.html")
library(htmltools)
library(webshot)
webshot::install_phantomjs()
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(widget_formattable, "./4-Correction/1000/Supplementary_table_3_CrypticMAPSs.png")
write.table(genes1000_corrected3, file="./4-Correction/1000/Genes_drivers3.csv", col.names=T, quote=F, row.names=F, sep = "\t")