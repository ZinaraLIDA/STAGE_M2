#!/bin/bash

#SBATCH --job-name=NeoRP22
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs_2

cd /home/lidamaha/NeoPredPipe/

python NeoRecoPo.py --neopred_in=/home/lidamaha/1-StageM2/8-Neoantigenes/NeoPredPipe/Sortie/b2_2/TestRun.neoantigens.txt --neoreco_out=/home/lidamaha/1-StageM2/8-Neoantigenes/NeoPredPipe/Sortie_NeoRecoPo/b2_2/ --fastas=/home/lidamaha/1-StageM2/8-Neoantigenes/NeoPredPipe/Sortie/b2_2/fastaFiles/
