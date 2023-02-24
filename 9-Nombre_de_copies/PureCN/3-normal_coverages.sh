#!/bin/bash

#SBATCH --job-name=norm_cov
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

ls -a $OUT/b1/*_loess.txt.gz | cat > $DIR/3-normal_coverages/b1/normal_coverages_b1.list
ls -a $OUT/b2/*_loess.txt.gz | cat > $DIR/3-normal_coverages/b2/normal_coverages_b2.list
