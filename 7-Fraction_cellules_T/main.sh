#!/bin/bash

#SBATCH --job-name=tcell
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs

/home/lidamaha/miniconda3/envs/ngs/bin/Rscript /home/lidamaha/1-StageM2/7-Fraction_cellules_T/TcellExTRECT.R
