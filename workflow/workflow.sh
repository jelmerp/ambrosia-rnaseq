## Scripts
scr_mqc=mcic-scripts/misc/fastqc.sh
scr_trimg=mcic-scripts/misc/trimgalore.sh
scr_idx=mcic-scripts/rnaseq/star_index.sh
scr_aln=mcic-scripts/rnaseq/star_align.sh

## Data and results
dir_scratch=/fs/scratch/PAS0471/jelmer/assist/2021-09_nisha

#subvar="/subset"
subvar=""

dir_fq=data/fastq"$subvar"
dir_fqc=results/fastqc"$subvar"
dir_multiqc=results/multiqc"$subvar"
dir_trim=results/trimgalore"$subvar"
dir_map2fun=results/map2fungal"$subvar"
dir_map2ambrox=results/map2ambrox"$subvar"
dir_rcorr=results/rcorrector"$subvar"
dir_readfilt=results/readfilter"$subvar"
dir_trinity_scratch="$dir_scratch"/results/trinity"$subvar"
dir_trinity=results/trinity"$subvar"
dir_cdhit=results/cdhit
dir_busco=results/busco
dir_rnaquast=results/rnaquast
dir_detonate=results/detonate

readstats=results/qc_misc/"$subvar"readstats.txt

## Reference data
dir_funref=results/fungal_ref
fun_refmix="$dir_funref"/fun_refmix.fna
fun_refmix_idx="$dir_funref"/star_idx

ref_ambrox=results/refdata/ambro_xylebori/GCA_002778035.1_ASM277803v1_genomic.fna.gz
ref_ambrox_idx=results/refdata/ambro_xylebori/star_idx

## Settings
idx_size_fun=14      # STAR index size for combined fungal ref "genome" -- grep -v "^>" $fungal_refmix | wc -c # 748503311 => log2(748503311)/2 - 1) = 13.74 = 14 
idx_size_ambro=11    # STAR index size for Ambrosiella ref genome -- zgrep -v "^>" $ref_ambro | wc -c # 27541947 => log2(27541947)/2 - 1) = 11.35 = 11
minqual=5            # Min. qual for TrimGalore trimming - follows https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
minlen=36            # Min. read length for TrimGalore trimming - follows https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.htmlq
mapmax_fun=100       # Max. nr of multimapping for fungal (pre-)mapping


# PREPROCESSING ----------------------------------------------------------------
## Subset FASTQ files for testing
bash mcic-scripts/misc/subsample_fq_dir.sh -i data/fastq -o data/fastq/subset -n 100000

## Run FastQC
for fastq_file in "$dir_fq"/*fastq.gz; do
    echo -e "\nFASTQ file: $fastq_file"
    sbatch "$scr_mqc" "$fastq_file" "$dir_fqc"/raw
done

## Run MultiQC
sbatch "$scr_mqc" "$dir_fqc" "$dir_multiqc"

## Error correction with rcorrector
sbatch scripts/rcorrector.sh "$dir_fq" "$dir_rcorr"

## Remove reads with unfixable errors, according to rcorrector
for R1 in "$dir_rcorr"/*R1*.cor.fq.gz; do
    sbatch scripts/readfilter.sh "$R1" "$dir_readfilt"
done

## Run TrimGalore
for R1 in "$dir_readfilt"/*R1*fastq.gz; do
    sbatch "$scr_trimg" -i "$R1" -o "$dir_trim"/rcorr \
        -O "$dir_fqc"/trimmed_rcorr -q "$minqual" -l "$minlen"
done
for R1 in "$dir_fq"/*R1*fastq.gz; do
    sbatch "$scr_trimg" -i "$R1" -o "$dir_trim"/nocorr \
        -O "$dir_fqc"/trimmed_nocorr -q "$minqual" -l "$minlen"
done

# MAP TO FUNGAL GENOMES --------------------------------------------------------
## Download fungal seqs and create indices
sbatch scripts/make_fungal_ref.sh "$dir_funref"/raw "$fun_refmix"
sbatch "$scr_idx" -i "$fun_refmix" -o "$fun_refmix_idx" -s "$idx_size_fun"
sbatch "$scr_idx" -i "$ref_ambrox" -o "$ref_ambrox_idx" -s "$idx_size_ambro"

## Map FASTQs to fungal refs
for R1 in "$dir_trim"/rcorr/*_R1_*fastq.gz; do
    sbatch "$scr_aln" -i "$R1" -o "$dir_map2fun"/rcorr -r "$fun_refmix_idx"  -m "$mapmax_fun"
done

for R1 in "$dir_trim"/notcorr/*_R1_*fastq.gz; do #TODO RUN THIS
    sbatch "$scr_aln" -i "$R1" -o "$dir_map2fun"/notcorr -r "$fun_refmix_idx"  -m "$mapmax_fun"
done

## Map FASTQs to Ambrosiella xylebori only
for fq_R1 in "$dir_trim"/*_R1_*fastq.gz; do
    sbatch "$scr_aln" -i "$fq_R1" -r "$ref_ambrox_idx" -o "$dir_map2ambrox" -m 10
done

## Create a table with stats on number of reads per step
bash scripts/readstats.sh "$dir_trim" "$dir_map2fun" "$dir_map2ambrox" "$readstats"


# ASSEMBLE BEETLE READS --------------------------------------------------------
#> https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html

## Assemble with Trinity
sbatch scripts/trinity.sh "$dir_map2fun"/rcorr/unmapped "$dir_trinity_scratch"
cp "$dir_trinity_scratch"/Trinity.fasta "$dir_trinity"

## Redundant transcript removal
sbatch scripts/cd-hit.sh $dir_trinity/Trinity.fasta $dir_cdhit/trinity_cl0.95.fa

## Run Trinity-stats - https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
TrinityStats.pl $dir_trinity/Trinity.fasta > "$dir_trinity"/Trinity_assembly.metrics

## rnaQUAST
sbatch scripts/rnaquast.sh $dir_cdhit/trinity_cl0.95.fa "$dir_rnaquast"

## Run BUSCO
busco_db=insecta_odb10
sbatch scripts/busco.sh "$dir_cdhit"/trinity_cl0.95.fa "$dir_busco" "$busco_db"
mv busco_*log "$dir_busco"; mv busco_downloads "$dir_busco"
        # --------------------------------------------------
        # |Results from dataset insecta_odb10               |
        # --------------------------------------------------
        # |C:93.6%[S:57.8%,D:35.8%],F:2.6%,M:3.8%,n:1367    |
        # |1279   Complete BUSCOs (C)                       |
        # |790    Complete and single-copy BUSCOs (S)       |
        # |489    Complete and duplicated BUSCOs (D)        |
        # |35     Fragmented BUSCOs (F)                     |
        # |53     Missing BUSCOs (M)                        |
        # |1367   Total BUSCO groups searched               |

## Run Detonate
sbatch scripts/detonate.sh "$dir_cdhit"/trinity_cl0.95.fa "$dir_trim"/rcorr "$dir_detonate"

## Quantify read support by mapping with Bowtie - https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
# bowtie2-build -–threads 4 $1 $assembly_prefix
# bowtie2 -p 10 -q --no-unal -k 20 -x $1 -1 $2 -2 $3  2>align_stats.txt| ${SINGULARITY_EXEC} samtools view -@10 -Sb -o bowtie2.bam

#! https://www.protocols.io/view/de-novo-transcriptome-assembly-workflow-ghebt3e?step=15

#TODO - Filter using expression level https://www.biorxiv.org/content/10.1101/035642v2.full ?
#> In general, for low coverage datasets (less than 20 million reads), filtering based on expression,
#> using TPM=1 as a threshold performs well, with Transrate filtering being too aggressive.
#> With higher coverage data (more than 60 million reads) Transrate filtering may be optimal, as may gene expression filtering using a threshold of TPM=0.5.

## Another round of fungal sequence removal?
#> See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6463014/, which used competitive BlastP


# BEETLE ASSEMBLY ANNOTATION ---------------------------------------------------


# AMBROSIELLA ------------------------------------------------------------------
## Assemble
#> https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity
#> If your transcriptome RNA-seq data are derived from a gene-dense compact genome, such as from fungal genomes,
#> where transcripts may often overlap in UTR regions, you can minimize fusion transcripts by leveraging the '--jaccard_clip' option if you have paired reads. 

