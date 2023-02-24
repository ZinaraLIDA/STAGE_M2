#!/bin/bash

#SBATCH --job-name=cov_b1
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs

export PURECN="/home/lidamaha/miniconda3/envs/ngs/lib/R/library/PureCN/extdata"
export OUT_REF="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/1-references"
export OUT="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/2-coverage"
export DIR="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline"

Rscript $PURECN/Coverage.R \
	--out-dir $OUT/b1 \
	--bam $DIR/batch1.list \
	--intervals $OUT_REF/intervals.txt \
	--cores 4
