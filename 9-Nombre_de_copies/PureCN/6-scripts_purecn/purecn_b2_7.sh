#!/bin/bash

#SBATCH --job-name=prcnb2_7
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs

export PURECN="/home/lidamaha/miniconda3/envs/ngs/lib/R/library/PureCN/extdata"
export OUT_REF="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/5-normaldb/references_b2_7"
export OUT="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/7-sorties_purecn"
export DIR="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline"
export SAMPLEID="D210328"

mkdir $OUT/$SAMPLEID
Rscript $PURECN/PureCN.R \
	--out $OUT/$SAMPLEID \
	--tumor $DIR/2-coverage/b2/${SAMPLEID}_S4_coverage_loess.txt.gz \
	--sampleid $SAMPLEID \
	--vcf /home/lidamaha/1-StageM2/1-Appel_de_variant/6-Avec_Mutect2/batch2/fichiersGZ/${SAMPLEID}.vcf \
	--normaldb $OUT_REF/normalDB_agilent_v8_hg38.rds \
	--intervals $DIR/1-references/intervals.txt \
	--genome hg38
