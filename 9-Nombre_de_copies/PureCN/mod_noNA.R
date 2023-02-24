setwd("~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/PureCN")
library(vcfR)

ids <- c("D181208", "D181210", "D181211", "D181212", "D181213", "D181215", "D202061", "D202062", "D202063",
         "D202064", "D202065", "D202066", "D210285", "D210288", "D210289", "D210290", "D210294", "D210295",
         "D210296", "D210297", "D210327", "D210328", "D210332", "D210334", "D210335", "D210338", "D210339")

for (id in ids) {
  data_ <- read.vcfR(sprintf("./vcf_filtrE_MAJ/%s.vcf", id))
  meta <- data_@meta
  fix <- data.frame(data_@fix)
  gt <- data.frame(data_@gt)
  
  fix$FILTER <- "."
  fix$ID <- "."
  
  data_new <- cbind(fix, gt)
  
  write.table(data_new, file=sprintf("./vcf_filtrE_MAJ_noNA/%s.vcf", id),
              col.names=F, quote=F, row.names=F, sep = "\t")
}

