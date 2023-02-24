#!/bin/bash

#SBATCH --job-name=normb2_4
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs

export PURECN="/home/lidamaha/miniconda3/envs/ngs/lib/R/library/PureCN/extdata"
export OUT_REF="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline/5-normaldb/references_b2_4"
export DIR="/home/lidamaha/1-StageM2/9-Nombre_de_copies/PureCN/Pipeline"

Rscript $PURECN/NormalDB.R \
	--out-dir $OUT_REF \
	--coverage-files $DIR/3-normal_coverages/b2/normal_coverages_b2_4.list \
	--genome hg38 \
	--assay agilent_v8
