#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --out=slurm-make_fungal_refs-%j.out

## Software 
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate blast-env

## Bash strict mode
set -euo pipefail

## Command-line args
dir_refs=$1
refmix=$2

## Constants
FUNGUS_REF_TABLE_1=refdata/fungal/Blaz2019_TableS2_ed.txt
FUNGUS_REF_TABLE_2=refdata/fungal/PRJNA395605_AssemblyDetails.txt

## If necessary, create the output dir:
mkdir -p "$dir_refs"
refmix_dir=$(basename "$refmix")
mkdir -p "$refmix_dir"

## Report before starting the pogram:
echo
echo "## Starting script to download fungal reference sequences..."
date
echo "## Dir for raw downloaded fungal ref seqs: $dir_refs"
echo "## Output file with concatenated sequences: $refmix"
echo -e "------------------\n\n"

#! NOTE: JGI SEQUENCES WERE MANUALLY DOWNLOADED FIRST

## Download sequences in first file
while IFS=$'\t' read -r -u 9 species db accession; do
    if [ "$db" = "NCBI" ]; then
        echo -e "\n----------\n## Species: $species / Accession: $accession"
        ftp_dir=$(esearch -db assembly -query "$accession" | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank)
        fname=$(echo "$ftp_dir" | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/')
        url="$ftp_dir"/"$fname"
        echo "URL: $url"

        fname_local="$dir_refs"/"$fname"
        [[ -f "$fname_local" ]] && echo "## Local file already present, not downloading..."
        [[ ! -f "$fname_local" ]] && wget -O "$fname_local" "$url"
        [[ ! -f "$fname_local" ]] && echo "## ERROR: Local file not present after download attempt" #& exit 1
    fi
done 9<$FUNGUS_REF_TABLE_1

## Download sequences in second file
while IFS=$'\t' read -ru 9 accession discard_rest; do
    echo -e "\n----------\n## Accession: $accession"
    echo "$discard_rest" > /dev/null
    ftp_dir=$(esearch -db assembly -query "$accession" | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank)
    fname=$(echo "$ftp_dir" | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/')
    url="$ftp_dir"/"$fname"
    echo "URL: $url"

    fname_local="$dir_refs"/"$fname"
    [[ -f "$fname_local" ]] && echo "## Local file already present, not downloading..."
    [[ ! -f "$fname_local" ]] && wget -O "$fname_local" "$url"
    [[ ! -f "$fname_local" ]] && echo "## ERROR: Local file not present after download attempt" #& exit 1
done 9<$FUNGUS_REF_TABLE_2

## Concatenate references
zcat "$dir_refs"/*.f*a.gz > "$refmix"

## Report
echo -e "\n## Listing downloaded files:"
ls -lh "$dir_refs"
echo -e "\n## Listing reference mix:"
ls -lh "$refmix"

echo -e "\n## Done with script make_fungal_refs.sh"
date
