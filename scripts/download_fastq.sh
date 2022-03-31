#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --out=slurm-download-%j.out

## Bash strict mode
set -euo pipefail

## Report
echo "## Starting download..."
date

## Variables
outdir=$1

## Constants
URL=http://lennon.gnets.ncsu.edu/Patwa/Patwa_NVS112A_fastqfiles.tar
PASSWORD=cSArtpVF
USERNAME=patwa

## Create the output dir
mkdir -p "$outdir"

## Download
wget --user "$USERNAME" --password "$PASSWORD" -P "$outdir" "$URL"

## Extract
tar -xvf "$outdir"/"$(basename "$URL")" "$outdir"
mv "$outdir"/Patwa_fastqfiles/* "$outdir"
rmdir "$outdir"/Patwa_fastqfiles

## Report
echo "## Done"
date
