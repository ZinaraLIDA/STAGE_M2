rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/4-Matrice_atlas")

# Importation du jeu de données
###############################

data_ <- readRDS( file = "~/ZINARA/Gnomic/1-StageM2/0-Metadata/Normal_epithelial_buccal_mucosa.rds" )

library(Seurat)

# Dans au moins 5% des cellules
################################

### Méthode 1
expressed_genes <- CreateSeuratObject(counts = data_, min.cells = round((data_@Dim[2]*5)/100, 0))
expressed_genes <- FindVariableFeatures(expressed_genes, selection.method = "vst", nfeatures = Inf)
expressed_genes_list <- VariableFeatures(expressed_genes)
if (!dir.exists("1-Sortie")) { dir.create("1-Sortie") }
write.table(expressed_genes_list, file="./1-Sortie/genes_exprimes_dans_les_cellulles_epitheliales_de_la_muqueuse_orale.txt", col.names=F, quote=F, row.names=F)

### Méthode 2
data_df <- as.data.frame(data_)[1:(data_@Dim[1]), 1:(data_@Dim[2])]
data_df_N <- data_df
data_df_N[data_df_N >= 1] <- 1
data_df_N <- transform(data_df_N, Effectif = rowSums(data_df_N))
data_df_N <- transform(data_df_N, Pourcentage = round( (data_df_N$Effectif)/(data_@Dim[2])*100, 2 ))
expressed_genes2 <- subset( data_df_N, Pourcentage >= 5 )
expressed_genes_list2 <- row.names(expressed_genes)

### Vérification
length(intersect(expressed_genes_list, expressed_genes_list2)) == length(expressed_genes_list)

# Absent
########

# Dans moins de 1% des cellules
unexpressed_genes <- subset( data_df_N, as.numeric(Pourcentage) < 1 )
unexpressed_genes <- subset( data_df_N, as.numeric(Effectif) < round((data_@Dim[2]*1)/100, 0) )

unexpressed_genes_list <- row.names(unexpressed_genes)
write.table(unexpressed_genes_list, file="./1-Sortie/genes_non_exprimes_dans_les_cellulles_epitheliales_de_la_muqueuse_orale.txt", col.names=F, quote=F, row.names=F)
