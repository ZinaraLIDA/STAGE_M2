rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes/Analyse")

ids <- c("D181208", "D181210", "D181211", "D181212", "D181213", "D181215", "D202061", "D202062",
         "D202063", "D202064", "D202065", "D202066", "D210285", "D210288", "D210289", "D210290",
         "D210294", "D210295", "D210296", "D210297", "D210327", "D210328", "D210332", "D210334",
         "D210335", "D210338", "D210339")

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


df_all <- data.frame(Sample=character(), Total_mutation=numeric(), Mutation_giving_neoantigen=numeric(), 
                     p_Mutation_giving_neoantigen=numeric(), Total_neoantigen=numeric(), Neoantigen_WeakBinder=numeric(),
                     p_Neoantigen_WeakBinder=numeric(), Neoantigen_StrongBinder=numeric(), p_Neoantigen_StrongBinder=numeric(),
                     Novelty=numeric(), p_Novelty=numeric())

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
  
  neo_Ag_i <- subset(neo_Ag, Sample==idx)
  neo_Ag_summary_i <- subset(neo_Ag_summary, Sample==idx)
  
  id_data_i <- unique(data_i$ID)
  id_neo_Ag_i <- unique(neo_Ag_i$id)
  intersection_i <- intersect(id_data_i, id_neo_Ag_i)
  difference_i <- setdiff(id_data_i, id_neo_Ag_i)
  
  # Pourcentage de mutations donnant lieu à des néoantigènes
  p_i <- round(length(intersection_i) / length(id_data_i) * 100, 2)
  # Pourcentage de néoantigènes qui sont des ligants fort
  posx <- grep(idx, neo_Ag_summary_i$Sample)
  SB_i <- neo_Ag_summary_i$Total_SB[posx]
  Total_i <- neo_Ag_summary_i$Total[posx]
  p_i_i <- round( SB_i / Total_i * 100, 2)
  # Pourcentage de néoantigènes qui sont des ligants faible
  WB_i <- neo_Ag_summary_i$Total_WB[posx]
  p_i_i_i <- round( WB_i / Total_i * 100, 2)
  # Pourcentage de nouveaux néoantigènes
  new <- nrow(subset(neo_Ag_i, Novelty==1))
  p_i_i_i_i <- round( new / Total_i * 100, 2)
  
  df_all[dim(df_all)[1]+1,] <- c(idx, length(id_data_i), length(intersection_i), p_i, Total_i, WB_i, p_i_i_i, SB_i, p_i_i, new, p_i_i_i_i)
}

