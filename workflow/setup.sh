#!/bin/bash

## Create Conda environments
mamba env create --force --file envs/rcorrector.yml
mamba env create --force --file envs/trinity.yml
mamba env create --force --file envs/busco.yml
mamba env create --force --file envs/cd-hit.yml
mamba env create --force --file envs/transrate.yml
mamba env create --force --file envs/detonate.yml 

## Bellepheron
cd software
git clone https://github.com/JesseKerkvliet/Bellerophon
conda env create --file Bellerophon/envs/Bellerophon.yml
cp Bellerophon/Bellerophon.py /users/PAS0471/jelmer/miniconda3/envs/bellerophon/bin/

## Transrate
https://github.com/blahah/transrate/issues/218

## Transfer files
rsync -avz jp:/fs/project/PAS0471/jelmer/assist/2021-09_nisha/results /home/jelmer/Dropbox/mcic/assist/2021-09_nisha
rsync -avz /home/jelmer/Dropbox/mcic/assist/2021-09_nisha/refdata/*gz jp:/fs/project/PAS0471/jelmer/assist/2021-09_nisha/refdata

## Create software environments
module load python/3.6-conda5.2
conda create -n trimgalore-env -y -c bioconda trim-galore=0.6.7
#conda create -p /users/PAS0471/jelmer/.conda/envs/trimgalore-env -c bioconda trim-galore=0.6.7

## Download the FASTQ files from the NCSU sequencing core website
dir_fq=data/fastq
URL=http://lennon.gnets.ncsu.edu/Patwa/Patwa_NVS112A_fastqfiles.tar
sbatch scripts/download_fastq.sh $URL "$dir_fq"

tar -xvf "$dir_fq"/"$(basename "$URL")" "$dir_fq"
mv "$dir_fq"/Patwa_fastqfiles/* "$dir_fq"
rmdir "$dir_fq"/Patwa_fastqfiles

## Download the Ambrosiella cleistominuta reference genome
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/139/545/GCA_017139545.1_ASM1713954v1/GCA_017139545.1_ASM1713954v1_genomic.fna.gz
wget "$URL" -P results/refdata/ambro_cleistominuta/