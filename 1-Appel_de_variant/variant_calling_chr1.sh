#!/bin/bash

#SBATCH --job-name=nschr1
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL


source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs

module load singularity

cd /home/lidamaha/1-StageM2/1-Appel_de_variant/1-Fichiers_bam_splitEs/1/

/home/lidamaha/nextflow run iarcbioinfo/needlestack \
	-with-singularity \
	./needlestack_latest.sif \
	--bed /home/martinep/data/Agilent_v6UTR_v8_intersect2_covered.bed \
	--input_bams ./chr1/ \
	--ref ./hg38-au.pri.fa \
	--output_vcf variants_chr1.vcf \
	-w ./
