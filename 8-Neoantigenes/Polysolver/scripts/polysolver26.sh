#!/bin/bash

#SBATCH --job-name=pls26
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com

source /home/lidamaha/miniconda3/etc/profile.d/conda.sh

conda activate polysolver

source /home/lidamaha/polysolver/scripts/config.sh
PERL5LIB=/home/lidamaha/miniconda3/envs/polysolver/lib/perl5/site_perl/5.22.0/ \
	/home/lidamaha/polysolver/scripts/shell_call_hla_type \
	/home/ferraria/work/dysect/data/seqruns/210607_A00317_0267_BHWYYHDMXX/aligned/bam/D210338_S12.bam \
	Unknown 1 hg38 STDFQ 0 \
	/home/lidamaha/1-StageM2/8-Neoantigenes/Polysolver/HLAtypes/D210338

