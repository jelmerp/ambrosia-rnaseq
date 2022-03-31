# DATA AND REF DATA ------------------------------------------------------------
## Download the FASTQ files from the NCSU sequencing core website
dir_fq=data/fastq
sbatch scripts/download_fastq.sh "$dir_fq"

## Download the Ambrosiella cleistominuta reference genome
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/139/545/GCA_017139545.1_ASM1713954v1/GCA_017139545.1_ASM1713954v1_genomic.fna.gz
wget "$URL" -P results/refdata/ambro_cleistominuta/


# SOFTWARE ---------------------------------------------------------------------
## Create Conda environments
mamba env create -f envs/rcorrector.yml
mamba env create -f envs/trinity.yml
mamba env create -f envs/busco.yml
mamba env create -f envs/cd-hit.yml
mamba env create -f envs/transrate.yml
mamba env create -f envs/detonate.yml 
mamba env create -f envs/trinonate.yml
mamba env create -f envs/kallisto.yml

## Install TRINOTATE dependencies into conda environment
cd software
wget https://services.healthtech.dtu.dk/download/b3dac17d-cc0b-4cf8-9ee7-92aafe1cb536/signalp-4.1g.Linux.tar.gz
tar -xzvf signalp-4.1g.Linux.tar.gz
cp signalp-4.1/signalp /users/PAS0471/jelmer/miniconda3/envs/trinotate-env/bin/

wget https://services.healthtech.dtu.dk/download/3c40df53-f9c2-47c2-80a3-798c06b5cf85/tmhmm-2.0c.Linux.tar.gz
tar -xzvf tmhmm-2.0c.Linux.tar.gz
cp tmhmm-2.0c/bin/* /users/PAS0471/jelmer/miniconda3/envs/trinotate-env/bin/

wget https://services.healthtech.dtu.dk/download/50fa2067-d5b4-409d-93d9-178b011680ed/rnammer-1.2.src.tar.gz
tar -xzvf rnammer-1.2.src.tar.gz
cp rnammer/*rnammer /users/PAS0471/jelmer/miniconda3/envs/trinotate-env/bin/

wget http://eddylab.org/software/hmmer/hmmer-2.3.tar.gz
tar -xzvf hmmer-2.3.tar.gz
cd hmmer-2.3
./configure
make
mv src/hmmsearch src/hmmsearch2

env_dir=/users/PAS0471/jelmer/miniconda3/envs/trinotate-env
mkdir "$env_dir"/bin/util/
cp "$env_dir"/bin/superScaffoldGenerator.pl "$env_dir"/bin/util/  # Trinotate will look for it in util
