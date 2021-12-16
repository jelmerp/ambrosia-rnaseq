#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=detonate
#SBATCH --output=slurm-detonate-%j.out

## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate detonate-env

## Bash strict mode
set -euo pipefail

## Command-line args
fa_in="$1"
indir="$2"
outdir="$3"

## Chck input
[[ "$#" != 3 ]] && echo "## ERROR: Please provide 3 arguments - you provided $#" && exit 1
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir $indir does not exist" && exit 1

## Process args
### Comma-delimited list of FASTQ files:
R1_list=$(echo "$indir"/*R1*fastq.gz | sed 's/ /,/g')
R2_list=$(echo "$indir"/*R2*fastq.gz | sed 's/ /,/g')

## Other paramaters
n_cores=$SLURM_CPUS_PER_TASK
readlen=150

## Report
echo
echo "## Starting script detonate.sh"
date
echo "## Args: $*"
echo
echo "## Input FASTA file:  $fa_in"
echo "## Output dir:        $outdir"
echo
echo "## List of R1 files:  $R1_list"
echo "## List of R2 files:  $R2_list"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN DETONATE -----------------------------------------------------------------
echo "## Starting Detonate run..."

echo -e "\n## Running rsem-eval-calculate-score..."
rsem-eval-calculate-score \
    --paired-end \
    "$R1_list" "$R2_list" \
    "$fa_in" \
    "$outdir" \
    "$readlen" \
    -p "$n_cores"

#TODO include transcript length reference
# --transcript-length-parameters rsem-eval/true_transcript_length_distribution/mouse.txt \
# http://deweylab.biostat.wisc.edu/detonate/vignette.html 


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script detonate.sh"
date


# INFO -------------------------------------------------------------------------
#? paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4298084/
#? website: http://deweylab.biostat.wisc.edu/detonate/vignette.html
#? conda: https://anaconda.org/bioconda/detonate

#> RSEM-EVAL is a reference-free evaluation method based on a novel probabilistic
#> model that depends only on an assembly and the RNA-Seq reads used for its construction.
#> Unlike N50, RSEM-EVAL combines multiple factors, including the compactness of an assembly
#> and the support of the assembly from the RNA-Seq data, into a single, statistically-principled evaluation score.
#> This score can be used to select a best assembler, optimize an assembler's parameters, and guide new assembler design as an objective function.
#> In addition, for each contig within an assembly, RSEM-EVAL provides a score that
#> assesses how well that contig is supported by the RNA-Seq data and can be used to filter unnecessary contigs.


