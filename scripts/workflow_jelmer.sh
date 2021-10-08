# working dir should be the project dir root
# rsync -avz jp:/fs/project/PAS0471/jelmer/assist/2021-09_nisha/results /home/jelmer/Dropbox/mcic/assist/2021-09_nisha

## Create software environments
module load python/3.6-conda5.2
conda create -n trimgalore-env -y -c bioconda trim-galore=0.6.7
#conda create -p /users/PAS0471/jelmer/.conda/envs/trimgalore-env -c bioconda trim-galore=0.6.7

## Set dirs
dir_fastq=data/fastq
dir_fastqc=results/fastqc
dir_multiqc=results/multiqc
dir_trim=results/trimgalore

## Download the FASTQ files
dir_fastq=data/fastq
URL=http://lennon.gnets.ncsu.edu/Patwa/Patwa_NVS112A_fastqfiles.tar
sbatch scripts/download_fastq.sh $URL $dir_fastq
URL=http://lennon.gnets.ncsu.edu/Patwa/EMySc1_S23_L002_R1_001.fastq.gz
sbatch scripts/download_fastq.sh $URL $dir_fastq

tar -xvf "$dir_fastq"/"$(basename "$URL")" "$dir_fastq"
mv "$dir_fastq"/Patwa_fastqfiles/* "$dir_fastq"
rmdir "$dir_fastq"/Patwa_fastqfiles

## Run FastQC
for fastq_file in "$dir_fastq"/*fastq.gz; do
    echo "file: $fastq_file"
    sbatch mcic-scripts/misc/fastqc.sh "$fastq_file" "$dir_fastqc"
done

## Run MultiQC
sbatch mcic-scripts/misc/multiqc.sh "$dir_fastqc" "$dir_multiqc"

## Run TrimGalore
for trimgalore_input in "$dir_fastq"/*fastq.gz; do
    echo -e "\nInput file: $trimgalore_input"
    sbatch scripts/trimgalore.sh "$trimgalore_input" "$dir_trim" "$dir_fastqc"
done