#!/bin/bash

#SBATCH --job-name=neop_b22
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs_2

cd /home/lidamaha/NeoPredPipe/

python NeoPredPipe.py -I /home/lidamaha/1-StageM2/8-Neoantigenes/vcf_filtrE/b2_2/ \
	-H /home/lidamaha/1-StageM2/8-Neoantigenes/Polysolver/HLAtypes/hlatypes_batch2_2.txt \
	-o /home/lidamaha/1-StageM2/8-Neoantigenes/NeoPredPipe/Sortie/b2_2/ \
	-n TestRun -c 0 -E 8 9 10 -d -m
