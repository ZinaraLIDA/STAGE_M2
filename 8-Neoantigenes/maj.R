library(vcfR)
setwd("~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes")


ids <- c("D181208", "D181210", "D181211", "D181212", "D181213", "D181215", "D202061", "D202062",
        "D202063", "D202064", "D202065", "D202066", "D210285", "D210288", "D210289", "D210290",
        "D210294", "D210295", "D210296", "D210297", "D210327", "D210328", "D210332", "D210334",
        "D210335", "D210338", "D210339")


for (id in ids) {
  data_ <- read.vcfR(sprintf("./vcf_filtrE/%s.vcf", id))
  meta <- data_@meta
  fix <- data.frame(data_@fix)
  gt <- data.frame(data_@gt)
  
  Y <- gt$FORMAT[1]
  mja_format <- function(x) {
    left <- substr(x, 3, nchar(x))
    all_ <- paste("GT", ":AD", left, sep = "")
    return(all_)
  }
  mja_format(Y)
  
  Z <- gt$Sample[1]
  mja_sample <- function(x) {
    first <- substr(x, 1, 3)
    left <- substr(x, 5, nchar(x))
    RO <- unlist(strsplit(x, ":"))[4]
    AO <- unlist(strsplit(x, ":"))[5]
    AD <- paste(RO, AO, sep = ",")
    all_ <- paste(first, AD, left, sep = ":")
    return(all_)
  }
  mja_sample(Z)
  
  new_gt <- gt
  new_gt$FORMAT <- unlist(lapply( 1:nrow(new_gt), function(i) { mja_format(new_gt$FORMAT[i]) } ))
  new_gt$Sample <- unlist(lapply( 1:nrow(new_gt), function(i) { mja_sample(new_gt$Sample[i]) } ))
  
  new_data <- cbind(fix, new_gt)
  
  write.table(new_data, file=sprintf("./vcf_filtrE_MAJ/%s.vcf", id), col.names=T, quote=F, row.names=F, sep = "\t")
}
