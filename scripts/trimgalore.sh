#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-trimmomatic-%j.out

## Load software
module load python/3.6-conda5.2
source activate trimgalore-env

## Bash strict settings
set -euo pipefail

## Report
echo "## Starting script trimgalore.sh"
date

## Check number of command-line args
if [ "$#" -ne 3 ]; then
    echo "Usage: trimgalore.sh <input-file.fastq[.gz]> <output-dir-trim> <output-dir-fastqc>"
    echo "Exiting."
    exit 1
fi

## Command-line args
input="$1"
outdir_trim="$2"
outdir_fastqc="$3"

## Process variables and args
mkdir -p "$outdir_trim"
mkdir -p "$outdir_fastqc"

n_threads="$SLURM_CPUS_ON_NODE"

## Report
echo "## Input file:                    $input"
echo "## Output dir - trimmed FASTQs:   $outdir_trim"
echo "## Output dir - FastQC:           $outdir_fastqc"
echo -e "---------------------------\n\n"

## Run Trim-Galore
trim_galore \
    --quality 0 --length 20 \
    --gzip -j "$n_threads" \
    --output_dir "$outdir_trim" \
    --fastqc --fastqc_args "-t $n_threads --outdir $outdir_fastqc" \
    "$input"

## Move FASTQ file
file_id=$(basename "$input" .fastq.gz)
mv "$outdir_trim"/"$file_id"_trimmed.fq.gz "$outdir_trim"/"$file_id"_trimmed.fastq.gz

## Report
echo -e "\nDone with script trimgalore.sh"
date