## Create software environments
module load python/3.6-conda5.2
conda create -n trimgalore-env -y -c bioconda trim-galore=0.6.7

## Set dirs
dir_fastq=data/fastq
dir_fastqc=results/fastqc
dir_multiqc=results/multiqc
dir_trim=results/trimgalore

## Download the FASTQ files
sbatch scripts/download_fastq.sh

tar -xvf "$dir_fastq"/Patwa_NVS112A_fastqfiles.tar
mv Patwa_fastqfiles/* data/fastq/
rmdir Patwa_fastqfiles

## Run FastQC
for fastq_file in data/fastq/*fastq.gz; do
    echo "The FASTQ file is: $fastq_file"
    sbatch scripts/fastqc.sh $fastq_file results/fastqc
done
