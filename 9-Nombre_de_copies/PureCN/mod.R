
library(vcfR)

data_db <- read.vcfR("./dbsnp3/00-All.vcf.gz")

meta_db <- data_db@meta
fix_db <- data.frame(data_db@fix)
gt_db <- data.frame(data_db@gt)

data_ <- read.vcfR("./vcf_filtrE_MAJ/D181208.vcf")
meta <- data_@meta
fix <- data.frame(data_@fix)
gt <- data.frame(data_@gt)

fix$FILTER <- "."
fix_db$FILTER <- "."
fix$ID <- "."
fix_db$ID <- "."

data_new <- cbind(fix, gt)
data_new["id"] <- paste(data_new$CHROM, data_new$POS, data_new$REF, data_new$ALT, sep = "")

fix_db$CHROM <- paste("chr", fix_db$CHROM, sep = "")
fix_db["id"] <- paste(fix_db$CHROM, fix_db$POS, fix_db$REF, fix_db$ALT, sep = "")

#union_ <- fix_db[fix_db$id %in% data_new$id, ]
data_new$ID <- unlist(lapply(1:nrow(data_new), function(i){fix_db$ID[grep(data_new$id[i], fix_db$id)]}))
data_new <- data_new[,1:10]







write.table(data_new, file="./vcf_filtrE_MAJ_/D181208.vcf",
            col.names=F, quote=F, row.names=F, sep = "\t")