#!/bin/bash

#SBATCH --job-name=an_12
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zinaralidamahasolo@gmail.com


source /home/lidamaha/miniconda3/etc/profile.d/conda.sh
conda activate ngs

cd /home/lidamaha/ensembl-vep/

/home/lidamaha/ensembl-vep/vep \
	-i /home/lidamaha/1-StageM2/2-Filtre/1-Fichiers_filtrEs/echantillon_12-filtrE.txt \
	-o /home/lidamaha/1-StageM2/3-Annotation/2-Fichiers_annotEs/echantillon_12-annotE.txt \
	--fasta /home/lidamaha/1-StageM2/1-Appel_de_variant/1-Fichiers_bam_splitEs/1/hg38-au.pri.fa \
	--cache --species homo_sapiens -a GRCh38 --cache_version 105 --sift b --offline \
	--force_overwrite \
	--custom /home/martinep/.vep/customs/clinvar_20210102.vcf.gz,clinvar,vcf,overlap,0,CLNSIG,CLNSIGCONF,CLNDISDB,CLNREVSTAT \
	--custom /home/martinep/.vep/customs/gnomad.genomes.r3.0.pass.af-only.hg38.vcf.gz,gnomad,vcf,overlap,0,AC,AF,nhomalt \
	--variant_class 1 --polyphen b --gene_phenotype --regulatory --phased --total_length --numbers --vcf_info_field CSQ \
	--hgvs --hgvsg --shift_hgvs 1 --protein --symbol --ccds --uniprot --tsl --appris --canonical --biotype --domains \
	--xref_refseq --check_existing --max_af --pubmed --gencode_basic --exclude_predicted --per_gene \
	--pick_order canonical,appris,tsl,biotype,ccds,rank,length

