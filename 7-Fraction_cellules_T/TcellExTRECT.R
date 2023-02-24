rm(list = ls())
setwd("/home/lidamaha/1-StageM2/7-Fraction_cellules_T/")

# Définition des chemins et calcul des valeurs de couverture
library(TcellExTRECT)

patients <- read.table( file = "/home/lidamaha/1-StageM2/0-Metadata/patients.txt", header = TRUE, stringsAsFactors = FALSE )
library(comprehenr)
batch1_id <- to_vec( for (i in 1:6) for (j in 1:length( unlist( strsplit( patients$IDs[i], "," ) ) )) unlist( strsplit( patients$IDs[i], "," ) )[j] )
batch2_id <- to_vec( for (i in 7:15) for (j in 1:length( unlist( strsplit( patients$IDs[i], "," ) ) )) unlist( strsplit( patients$IDs[i], "," ) )[j] )
batch2 <- list.files("/home/ferraria/work/dysect/data/seqruns/210607_A00317_0267_BHWYYHDMXX/aligned/bam/")
library(stringr)
batch2 <- batch2[str_detect(string = batch2, pattern = ".bam$")]

TCRA_out_all <- data.frame( sample = character(), TCRA.tcell.fraction = numeric(), TCRA.tcell.fraction.lwr = numeric(), TCRA.tcell.fraction.upr = numeric(),
                            qcFit = numeric() )
idx = 1
for (id in batch1_id) {
  
  bam_file_path <- sprintf( '/home/martinep/data/dysect/%s_merged_markdup.bam', id )
  
  bed.file <- '/home/lidamaha/1-StageM2/1-Appel_de_variant/1-Fichiers_bam_splitEs/1/Agilent_v6UTR_v8_intersect2_covered.bed'
  data("tcra_seg_hg38")
  cov_example <- getCovFromBam(bamPath = bam_file_path,
                               outPath = './TCRA_files/',
                               vdj.seg = tcra_seg_hg38)
  cov_df <- loadCov(cov_example)
  
  # Exécution
  TCRA_exons_hg38_custom <- createExonDFBed(bed.file, 'hg38')
  TCRA.out <- runTcellExTRECT(cov_df, TCRA_exons_hg38_custom, tcra_seg_hg38, 'hg38')
  
  NewRow <- c(id, TCRA.out[[2]], TCRA.out[[3]], TCRA.out[[4]], TCRA.out[[5]])
  TCRA_out_all[nrow(TCRA_out_all)+1,] <- NewRow
  
  # Vérification graphique s'il n'y a pas d'exons ou de reads qui sont aberrants et peuvent interférer dans le calcul des fractions des cellules T
  png(file = sprintf("./images/%s-TCellExTRECT_log_ratio_%s.png", idx, id), width = 900, height = 800)
  plotTcellExTRECT(cov_df, TCRA_exons_hg38_custom, tcra_seg_hg38,'hg38', sample_name = 'TEST')
  dev.off()
  
  idx = idx+1
  
}
print(median(as.numeric(TCRA_out_all[substr(TCRA_out_all$sample,1,3)=="D18",]$TCRA.tcell.fraction)))
print(median(as.numeric(TCRA_out_all[substr(TCRA_out_all$sample,1,3)=="D20",]$TCRA.tcell.fraction)))
print(median(as.numeric(TCRA_out_all$TCRA.tcell.fraction)))

#####
for (i in 1:length(batch2)) {

  bam_file_path <- sprintf( '/home/ferraria/work/dysect/data/seqruns/210607_A00317_0267_BHWYYHDMXX/aligned/bam/%s', batch2[i] )

  bed.file <- '/home/lidamaha/1-StageM2/1-Appel_de_variant/1-Fichiers_bam_splitEs/1/Agilent_v6UTR_v8_intersect2_covered.bed'
  data("tcra_seg_hg38")
  cov_example <- getCovFromBam(bamPath = bam_file_path,
                               outPath = './TCRA_files/',
                               vdj.seg = tcra_seg_hg38)
  cov_df <- loadCov(cov_example)

  # Exécution
  TCRA_exons_hg38_custom <- createExonDFBed(bed.file, 'hg38')
  TCRA.out <- runTcellExTRECT(cov_df, TCRA_exons_hg38_custom, tcra_seg_hg38, 'hg38')

  NewRow <- c(batch2_id[i], TCRA.out[[2]], TCRA.out[[3]], TCRA.out[[4]], TCRA.out[[5]])
  TCRA_out_all[nrow(TCRA_out_all)+1,] <- NewRow

  # Vérification graphique s'il n'y a pas d'exons ou de reads qui sont aberrants et peuvent interférer dans le calcul des fractions des cellules T
  png(file = sprintf("./images/%s-TCellExTRECT_log_ratio_%s.png", idx, batch2_id[i]), width = 900, height = 800)
  plotTcellExTRECT(cov_df, TCRA_exons_hg38_custom, tcra_seg_hg38,'hg38', sample_name = 'TEST')
  dev.off()

  idx = idx+1

}
print(median(as.numeric(TCRA_out_all[substr(TCRA_out_all$sample,1,3)=="D21",]$TCRA.tcell.fraction)))

write.table(TCRA_out_all, file="./TCRA_out_all.txt", col.names=T, quote=F, row.names=F, sep = "\t")
