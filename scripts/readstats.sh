#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-readstats-%j.out

## Bash strict settings
set -euo pipefail

## Command-line args
dir_trim=$1
dir_map2fun=$2
dir_map2ambro=$3
readstats_file=$4

## Report
echo
echo "## Starting script readstats.sh"
date
echo "## Dir with trimmed FASTQs: $dir_trim"
echo "## Dir with BAM files mapped to all fungi: $dir_map2fun"
echo "## Dir with BAM files mapped to ambrosiella only: $dir_map2ambro"
echo "## Output file: $readstats_file"
echo -e "--------------------------\n"

## Create file header with column names
echo sample_id nread_raw nbase_raw nread_trim pbase_trim nread_map2fun nread_map2fun_unq nread_map2fun_multi pread_map2fun \
     nread_map2ambro nread_map2ambro_unq nread_map2ambro_multi pread_map2ambro |
     sed -E 's/ /\t/g' > "$readstats_file"
#echo "sample_id nread_raw nbase_raw nread_trim pbase_trim nread_map2fun pread_map2fun nread_map2ambro pread_map2ambro" | sed -E 's/ /\t/g' > "$readstats_file"

## Loop through samples
for fq in data/fastq/*R1*fastq.gz; do
    sam=$(basename "$fq" | sed 's/_.*//')

    nread_raw=$(grep "Total reads processed" "$dir_trim"/"$sam"*R1*trimming_report.txt | sed 's/[^0-9]//g')
    nbase_raw=$(grep "Total basepairs processed" "$dir_trim"/"$sam"*R1*trimming_report.txt | sed 's/[^0-9]//g')
    nread_trim=$(grep "Reads written" "$dir_trim"/"$sam"*R1*trimming_report.txt | sed -E 's/.*: +([0-9,]+) \(.*/\1/' | sed 's/,//g')
    nbase_trim=$(grep "Total written (filtered)" "$dir_trim"/"$sam"*R1*trimming_report.txt | sed -E 's/.*: +([0-9,]+) bp.*/\1/' | sed 's/,//g')
    pbase_trim=$(python -c "print(round(($nbase_raw - $nbase_trim) / $nbase_raw, 3))")
    
    nread_map2fun_unq=$(grep "Uniquely mapped reads number" "$dir_map2fun"/star_logs/"$sam"*Log.final.out | sed 's/[^0-9]//g')
    nread_map2fun_multi=$(grep "Number of reads mapped to multiple loci" "$dir_map2fun"/star_logs/"$sam"*Log.final.out | sed 's/[^0-9]//g')
    nread_map2fun=$(( nread_map2fun_unq + nread_map2fun_multi ))
    pread_map2fun=$(python -c "print(round($nread_map2fun / $nread_trim, 3))")

    nread_map2ambro_unq=$(grep "Uniquely mapped reads number" "$dir_map2ambro"/star_logs/"$sam"*Log.final.out | sed 's/[^0-9]//g')
    nread_map2ambro_multi=$(grep "Number of reads mapped to multiple loci" "$dir_map2ambro"/star_logs/"$sam"*Log.final.out | sed 's/[^0-9]//g')
    nread_map2ambro=$(( nread_map2ambro_unq + nread_map2ambro_multi ))
    pread_map2ambro=$(python -c "print(round($nread_map2ambro / $nread_trim, 3))")
    
    #echo "$sam $nread_raw" "$nbase_raw" "$nread_trim" "$pbase_trim" "$nread_map2fun" "$pread_map2fun" "$nread_map2ambro" "$pread_map2ambro" | sed -E 's/ /\t/g' >> "$readstats_file"
    echo "$sam $nread_raw" "$nbase_raw" "$nread_trim" "$pbase_trim" \
         "$nread_map2fun" "$nread_map2fun_unq" "$nread_map2fun_multi" "$pread_map2fun" \
         "$nread_map2ambro" "$nread_map2ambro_unq" "$nread_map2ambro_multi" "$pread_map2ambro" |
         sed -E 's/ /\t/g' >> "$readstats_file"

done

## Report
echo "## Showing output file:"
column -t "$readstats_file"

echo -e "\n## Done with scripts readstats.sh"
date
