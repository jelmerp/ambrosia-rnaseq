#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --output=slurm-fastqc-%j.out

set -euo pipefail

## Load the FastQC module
module load fastqc

## Process command-line arguments
fastq_file=$1
output_dir=$2

## Report
echo "Starting fastqc.sh script..."
date
echo "The FASTQ file is: $fastq_file"
echo "The output dir is: $output_dir"

## Create output dir if needed
mkdir -p $output_dir

## Run FastQC
fastqc --outdir=$output_dir $fastq_file

## Report
echo "Done with the fastqc script."
date
