#!/bin/bash

#SBATCH --job-name=int
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs

export PURECN="/home/lidamaha/miniconda3/envs/ngs/lib/R/library/PureCN/extdata"
export OUT_REF="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/1-references"

Rscript $PURECN/IntervalFile.R \
	--in-file /home/lidamaha/1-StageM2/1-Appel_de_variant/1-Fichiers_bam_splitEs/1/Agilent_v6UTR_v8_intersect2_covered.bed \
	--fasta /home/lidamaha/1-StageM2/1-Appel_de_variant/1-Fichiers_bam_splitEs/1/hg38-au.pri.fa \
	--out-file $OUT_REF/intervals.txt \
	--off-target --genome hg38 \
	--export $OUT_REF/optimized_hg38.bed \
	--mappability /home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw
