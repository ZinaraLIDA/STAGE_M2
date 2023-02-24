#!/bin/bash

#SBATCH --job-name=prcnb1_4
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs

export PURECN="/home/lidamaha/miniconda3/envs/ngs/lib/R/library/PureCN/extdata"
export OUT_REF="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/5-normaldb/references_b1_4"
export OUT="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/7-sorties_purecn"
export DIR="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline"
export SAMPLEID="D181212"

mkdir $OUT/$SAMPLEID
Rscript $PURECN/PureCN.R \
	--out $OUT/$SAMPLEID \
	--tumor $DIR/2-coverage/b1/${SAMPLEID}_merged_markdup_coverage_loess.txt.gz \
	--sampleid $SAMPLEID \
	--vcf /home/lidamaha/1-StageM2/1-Appel_de_variant/6-Avec_Mutect2/batch1/fichiersGZ/${SAMPLEID}.vcf \
	--normaldb $OUT_REF/normalDB_agilent_v6_hg38.rds \
	--intervals $DIR/1-references/intervals.txt \
	--genome hg38
