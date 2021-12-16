#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --out=slurm-featurecounts-%j.out

# Software
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source ~/.bashrc
source activate subread-env

## Bash strict mode
set -euo pipefail

## Command-line args
indir=$1                     # Dir with input BAM file (every BAM file will be used)
outfile=$2                   # Output file with gene counts
gff=$3                       # Reference genome GFF file
t_option=${4-gene}           # featureCounts -t option (feature type)
g_option=${5-Name}           # featureCounts -g option (aggregration ID)

## Process args
outdir=$(dirname "$outfile")
mkdir -p "$outdir"

# Report
echo "## Starting script feature-counts.sh"
date 
echo "## Input dir with BAM files: $indir"
echo "## Output file with gene counts: $outfile"
echo "## Reference GFF file: $gff"
echo "## -t option (feature type): $t_option"
echo "## -g option (aggregration ID): $g_option"
echo
echo "## Number of BAM files: $(find "$outdir"/*bam | wc -l)"
echo -e "----------------------\n"

# Run featureCounts
featureCounts -s 2 -p -B -C \
    -T "$SLURM_CPUS_PER_TASK" \
    -t "$t_option" -g "$g_option" \
    -a "$gff" \
    -o "$outfile" \
    "$outdir"/*bam

## Report
echo -e "\n-------------------------"
echo "## Listing output file:"
ls -lh "$outfile"
echo -e "\n## Done with script feature-counts.sh"
date