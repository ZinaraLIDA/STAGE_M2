rm(list = ls())
setwd("~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/")

library(PureCN)

# 1.1-Target information
reference.file <- "~/ZINARA/Gnomic/1-StageM2/1-Appel_de_variant/1-Fichiers_bam_splitEs/1/hg38-au.pri.fa"
bed.file <- "~/ZINARA/Gnomic/1-StageM2/1-Appel_de_variant/1-Fichiers_bam_splitEs/1/Agilent_v6UTR_v8_intersect2_covered.bed"
intervals <- import(bed.file)
mappability.file <- "~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/PureCN/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw"
mappability <- import(mappability.file)
# Interval.R
interval.file <- preprocessIntervals(intervals, reference.file,
                                    mappability=mappability,
                                     output.file = "ex2_gc_file.txt")

# 1.2-Coverage data
bam.file <- "~/ZINARA/Gnomic/AVEC_NEEDLESTACK_2/Fraction_cellules_T/D181208_merged_markdup.bam"
interval.file <- "~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/ex2_gc_file.txt"
# Coverage.R
calculateBamCoverageByInterval(bam.file = bam.file,
                               interval.file = interval.file, output.file = "ex1_coverage.txt")

# 1.3-Library-specific coverage bias
normal.coverage.files <- "~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/..."
tumor.coverage.file <- "~/ZINARA/Gnomic/1-StageM2/9-Nombre_de_copies/ex1_coverage.txt"
vcf.file <- "~/ZINARA/Gnomic/1-StageM2/8-Neoantigenes/vcf_files_gz/..."

# 1.4-PON
# NormalDB.R
normalDB <- createNormalDatabase(normal.coverage.files)
saveRDS(normalDB, file="normalDB.rds")
normalDB <- readRDS("normalDB.rds")
pool <- calculateTangentNormal(tumor.coverage.file, normalDB)

# 1.7-Recommended run
# PureCN.R
ret <-runAbsoluteCN(normal.coverage.file=pool, # normal.coverage.file=normal.coverage.file,
                    tumor.coverage.file=tumor.coverage.file,
                    vcf.file=vcf.file,
                    genome="hg38", sampleid="Sample1",
                    interval.file=interval.file,
                    normalDB=normalDB,
                    args.filterVcf=list(),
                    args.filterIntervals = list(min.total.counts = 50),
                    post.optimize=FALSE, plot.cnv=FALSE, verbose=FALSE)








mappability.file <- system.file("extdata", "ex2_mappability.bigWig",
                                package = "PureCN", mustWork = TRUE)

mappability <- import(mappability.file)
