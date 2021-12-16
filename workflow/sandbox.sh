#!/bin/bash

## QC and filtering with TransRate
dir_transrate=results/transrate
sbatch -t 60 -c 1 --mem=4G scripts/transrate.sh $dir_cdhit/trinity_cl0.95.fa "$dir_trim"/rcorr "$dir_transrate"
#! FAILS! After fixing installation errors, now the SNAP aligner within TransRate keeps failing 

## Map FASTQs to Ambrosiella cleistominuta only
dir_trim=results/trimgalore
idx_size_ambro=11

ref_ambroc=results/refdata/ambro_cleistominuta/GCA_017139545.1_ASM1713954v1_genomic.fna.gz
ref_ambroc_idx=results/refdata/ambro_cleistominuta/star_idx
dir_map2ambroc=results/map2ambroc

sbatch mcic-scripts/rnaseq/star_index.sh -i "$ref_ambroc" -o "$ref_ambroc_idx" -s "$idx_size_ambro"

for fq_R1 in "$dir_trim"/*_R1_*fastq.gz; do
    sbatch mcic-scripts/rnaseq/star_align.sh -i "$fq_R1" -r "$ref_ambroc_idx" -o "$dir_map2ambroc" -m 100
done

#? > RESULT: ~10% lower number of reads map two cleistominuta vs xylebori
# (star-env) [jelmer@pitzer-login03 2021-09_nisha]$ grep "Uniquely mapped reads number" results/map2ambrox/star_logs/*Log.final.out
# results/map2ambrox/star_logs/EMySc1_S23_L002_Log.final.out:                   Uniquely mapped reads number |    3042417
# results/map2ambrox/star_logs/EMySc2_S7_L002_Log.final.out:                   Uniquely mapped reads number |     4118041
# results/map2ambrox/star_logs/EMySc3_S8_L002_Log.final.out:                   Uniquely mapped reads number |     3561314
# results/map2ambrox/star_logs/EMySc4_S9_L002_Log.final.out:                   Uniquely mapped reads number |     2661143
# results/map2ambrox/star_logs/FMy3_S12_L002_Log.final.out:                   Uniquely mapped reads number |      14342649
# results/map2ambrox/star_logs/FMy4_S13_L002_Log.final.out:                   Uniquely mapped reads number |      5264053
# results/map2ambrox/star_logs/FMySc2_S3_L002_Log.final.out:                   Uniquely mapped reads number |     16266698
# results/map2ambrox/star_logs/FMySc3_S4_L002_Log.final.out:                   Uniquely mapped reads number |     12505984

# (star-env) [jelmer@pitzer-login03 2021-09_nisha]$ grep "Uniquely mapped reads number" results/map2ambroc/star_logs/*Log.final.out
# results/map2ambroc/star_logs/EMySc1_S23_L002_Log.final.out:                   Uniquely mapped reads number |    2762826
# results/map2ambroc/star_logs/EMySc2_S7_L002_Log.final.out:                   Uniquely mapped reads number |     3805152
# results/map2ambroc/star_logs/EMySc3_S8_L002_Log.final.out:                   Uniquely mapped reads number |     3311184
# results/map2ambroc/star_logs/EMySc4_S9_L002_Log.final.out:                   Uniquely mapped reads number |     2437742
# results/map2ambroc/star_logs/FMy1_S10_L002_Log.final.out:                   Uniquely mapped reads number |      22464881
# results/map2ambroc/star_logs/FMy2_S11_L002_Log.final.out:                   Uniquely mapped reads number |      22051969
# results/map2ambroc/star_logs/FMy3_S12_L002_Log.final.out:                   Uniquely mapped reads number |      12880593
# results/map2ambroc/star_logs/FMy4_S13_L002_Log.final.out:                   Uniquely mapped reads number |      4814151
# results/map2ambroc/star_logs/FMySc1_S2_L002_Log.final.out:                   Uniquely mapped reads number |     15078761
# results/map2ambroc/star_logs/FMySc2_S3_L002_Log.final.out:                   Uniquely mapped reads number |     14611973
# results/map2ambroc/star_logs/FMySc3_S4_L002_Log.final.out:                   Uniquely mapped reads number |     11455483
# results/map2ambroc/star_logs/FMySc4_S5_L002_Log.final.out:                   Uniquely mapped reads number |     16878865