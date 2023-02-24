
rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/2-Filtre")
library(vcfR)

T1<-Sys.time()
###########################################################################################
###         Objectif 1 : VAF < 0.3, nbReads > 3, profondeur > 50, GT = "0/1"          #####
###########################################################################################

maxAF = 0.3
minReads = 3
minDP = 50
minQUALb1 = 100
minQUALb2 = 400


##############

patients <- read.table( file = "~/ZINARA/Gnomic/1-StageM2/0-Metadata/patients.txt", header = TRUE, stringsAsFactors = FALSE )
files <- list.files("../1-Appel_de_variant/4-Fichiers_variants_splitEs/")
ID <- unlist( lapply( 1:nrow(patients), function(i) { unlist( strsplit( patients$IDs[i], "," ) ) } ) )
fichiers_par_echantillon <- lapply( 1:length(ID), function(i) { files[ grep( ID[i], files) ] } )

# Pour batch1
#############

b1_ <- data.frame( "CHROM"=NA,"POS"=NA,"ID"=NA,"REF"=NA,"ALT"=NA,"QUAL"=NA,"FILTER"=NA,"INFO"=NA,"FORMAT"=NA,
                   "Details"=NA, "AF"=NA, "Sample"=NA, "Patient"=NA, stringsAsFactors = FALSE )
b1_ <- b1_[0,]
for (i in 1:12) {
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
      # On filtre DP AF  RO AO et QUAL, et on ne garde que les hétérozygotes (carrying 1 copy of each of the REF and ALT alleles)
      print( "     On filtre DP AF RO AO et GT" )
      data_filtrE_vcf <- data_[ which( (extract.gt(x = data_, element = "DP", as.numeric = TRUE) >= minDP ) &
                                         (extract.gt(x = data_, element = "AF", as.numeric = TRUE) <= maxAF ) &
                                         (extract.gt(x = data_, element = "RO", as.numeric = TRUE) >= minReads ) &
                                         (extract.gt(x = data_, element = "AO", as.numeric = TRUE) >= minReads ) &
                                         (extract.gt(x = data_, element = "GT", as.numeric = TRUE) == 0 ) ) ]
      print( "     On filtre QUAL" )
      data_filtrE_vcf <- data_filtrE_vcf[ which( getQUAL(data_filtrE_vcf) >= minQUALb1 ) ]
      # Vérification s'il reste des mutations après filtre
      print( "          Vérification s'il reste des mutations après filtre" )
      if ( length( data_filtrE_vcf@fix ) != 0 ) {
        # On remplit la dataframe pour l'échantillon i
        data_filtrE_df <- data.frame( data_filtrE_vcf@fix, data_filtrE_vcf@gt, 
                                      extract.gt(x = data_filtrE_vcf, element = "AF", as.numeric = TRUE), stringsAsFactors = FALSE )
        colnames(data_filtrE_df) <- c( "CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "Details", "AF" )
        # On crée une colonne Sample et une colonne Patient
        data_filtrE_df["Sample"] <- rep( ID[i], nrow(data_filtrE_df) )
        data_filtrE_df["Patient"] <- rep( patients$numero[grep(ID[i], patients$IDs)], nrow(data_filtrE_df) )
        dataframe_echantillon_i <- rbind(dataframe_echantillon_i, data_filtrE_df, stringsAsFactors=F)
        # On remplit la colonne ID avec CHRM+POS+REF+ALT+Sample+Patient
        print( "On remplit la colonne ID avec CHRM+POS+REF+ALT+Sample+Patient" )
        dataframe_echantillon_i$ID <- paste( dataframe_echantillon_i$CHROM, dataframe_echantillon_i$POS, dataframe_echantillon_i$REF,
                                             dataframe_echantillon_i$ALT, dataframe_echantillon_i$Sample, dataframe_echantillon_i$Patient, sep = "_" )
      }
    }
  }
  # On crée un fichier des données filtrées
  print( "On crée un fichier des données filtrées" )
  write.table(dataframe_echantillon_i[, 1:10], file=sprintf("./1-Fichiers_filtrEs/echantillon_%s-filtrE.txt", i), col.names=F,
              quote=F, row.names=F, sep = "\t")
  write.table(dataframe_echantillon_i, file=sprintf("./2-Fichiers_filtrEs_avec_toutes_les_colonnes/echantillon_%s-filtrE.txt", i),
              col.names=F, quote=F, row.names=F, sep = "\t")
  #######################
  b1_ <- rbind( b1_, dataframe_echantillon_i )
}
b1_ <- b1_[order(b1_$CHROM, as.numeric(b1_$POS)),]
write.table(b1_[, 1:10], file="./1-Fichiers_filtrEs/batch1-filtrE.txt", col.names=F,
            quote=F, row.names=F, sep = "\t")

# Pour batch2
#############

b2_ <- data.frame( "CHROM"=NA,"POS"=NA,"ID"=NA,"REF"=NA,"ALT"=NA,"QUAL"=NA,"FILTER"=NA,"INFO"=NA,"FORMAT"=NA,
                   "Details"=NA, "AF"=NA, "Sample"=NA, "Patient"=NA, stringsAsFactors = FALSE )
b2_ <- b2_[0,]
for (i in 13:27) {
  # On met dans un seul dataframe les mutations filtrées de l'échantillon i
  dataframe_echantillon_i <- data.frame( "CHROM"=NA,"POS"=NA,"ID"=NA,"REF"=NA,"ALT"=NA,"QUAL"=NA,"FILTER"=NA,"INFO"=NA,"FORMAT"=NA,
                                         "Details"=NA, "AF"=NA, "Sample"=NA, "Patient"=NA, stringsAsFactors = FALSE )
  dataframe_echantillon_i <- dataframe_echantillon_i[0,]
  for (fichier in fichiers_par_echantillon[[i]]) {
    # Importation du jeu de donnée
    print( sprintf( "--------------------------Importation du jeu de données %s --------------------------------------", fichier ) )
    data_ <- read.vcfR( file = sprintf( "../1-Appel_de_variant/4-Fichiers_variants_splitEs/%s", fichier ), verbose = FALSE )
    # Vérification si l'échantillon présente des mutations
    print( "     Vérification si l'échantillon présente des mutations" )
    if ( length( data_@fix ) != 0 ) {
      # On filtre DP AF  RO AO et QUAL, et on ne garde que les hétérozygotes (carrying 1 copy of each of the REF and ALT alleles)
      print( "     On filtre DP AF AO et GT" )
      data_filtrE_vcf <- data_[ which( (extract.gt(x = data_, element = "DP", as.numeric = TRUE) >= minDP ) &
                                         (extract.gt(x = data_, element = "AF", as.numeric = TRUE) <= maxAF ) &
                                         #(extract.gt(x = data_, element = "RO", as.numeric = TRUE) >= minReads ) &
                                         (extract.gt(x = data_, element = "AO", as.numeric = TRUE) >= minReads ) &
                                         (extract.gt(x = data_, element = "GT", as.numeric = TRUE) == 0 ) ) ]
      print( "     On filtre QUAL" )
      data_filtrE_vcf <- data_filtrE_vcf[ which( getQUAL(data_filtrE_vcf) >= minQUALb2 ) ]
      # Vérification s'il reste des mutations après filtre
      print( "          Vérification s'il reste des mutations après filtre" )
      if ( length( data_filtrE_vcf@fix ) != 0 ) {
        # On remplit la dataframe pour l'échantillon i
        data_filtrE_df <- data.frame( data_filtrE_vcf@fix, data_filtrE_vcf@gt, 
                                      extract.gt(x = data_filtrE_vcf, element = "AF", as.numeric = TRUE), stringsAsFactors = FALSE )
        colnames(data_filtrE_df) <- c( "CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "Details", "AF" )
        # On crée une colonne Sample et une colonne Patient
        data_filtrE_df["Sample"] <- rep( ID[i], nrow(data_filtrE_df) )
        data_filtrE_df["Patient"] <- rep( patients$numero[grep(ID[i], patients$IDs)], nrow(data_filtrE_df) )
        dataframe_echantillon_i <- rbind(dataframe_echantillon_i, data_filtrE_df, stringsAsFactors=F)
        # On remplit la colonne ID avec CHRM+POS+REF+ALT+Sample+Patient
        print( "On remplit la colonne ID avec CHRM+POS+REF+ALT+Sample+Patient" )
        dataframe_echantillon_i$ID <- paste( dataframe_echantillon_i$CHROM, dataframe_echantillon_i$POS, dataframe_echantillon_i$REF,
                                             dataframe_echantillon_i$ALT, dataframe_echantillon_i$Sample, dataframe_echantillon_i$Patient, sep = "_" )
      }
    }
  }
  # On crée un fichier des données filtrées
  print( "On crée un fichier des données filtrées" )
  write.table(dataframe_echantillon_i[, 1:10], file=sprintf("./1-Fichiers_filtrEs/echantillon_%s-filtrE.txt", i), col.names=F,
              quote=F, row.names=F, sep = "\t")
  write.table(dataframe_echantillon_i, file=sprintf("./2-Fichiers_filtrEs_avec_toutes_les_colonnes/echantillon_%s-filtrE.txt", i),
              col.names=F, quote=F, row.names=F, sep = "\t")
  #######################
  b2_ <- rbind( b2_, dataframe_echantillon_i )
}
b2_ <- b2_[order(b2_$CHROM, as.numeric(b2_$POS)),]
write.table(b2_[, 1:10], file="./1-Fichiers_filtrEs/batch2-filtrE.txt", col.names=F,
            quote=F, row.names=F, sep = "\t")

b_all <- rbind( b1_, b2_ )
write.table(b_all[, 1:10], file="./1-Fichiers_filtrEs/echantillon_fusionnE-filtrE.txt",
            col.names=F, quote=F, row.names=F, sep = "\t")
write.table(b_all, file="./2-Fichiers_filtrEs_avec_toutes_les_colonnes/echantillon_fusionnE-filtrE.txt",
            col.names=F, quote=F, row.names=F, sep = "\t")
###########################################################################

T2<-Sys.time()
Tdiff= difftime(T2, T1)
print(Tdiff)

quantile(as.numeric(b1_$QUAL))
quantile(as.numeric(b2_$QUAL))