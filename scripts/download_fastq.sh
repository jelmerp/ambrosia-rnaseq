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
URL=$1
outdir=$2

## Create the output dir
mkdir -p "$outdir"

## Download
### Using the username and password sent to us by the NCSU staff
wget --user patwa --password cSArtpVF -P "$outdir" "$URL"

## Report
echo "## Done"
date
