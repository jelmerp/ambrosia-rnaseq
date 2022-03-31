## Scripts
scr_mqc=mcic-scripts/qc/fastqc.sh
scr_trimg=mcic-scripts/trim/trimgalore.sh
scr_idx=mcic-scripts/rnaseq/star_index.sh
scr_aln=mcic-scripts/rnaseq/star_align.sh
bin=mcic-scripts/trans-assembly

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
dir_kallisto=results/kallisto/beetle
dir_trinotate=$PWD/results/trinotate

assembly_id=trinity_cl0.95
assembly=$PWD/"$dir_cdhit"/trinity_cl0.95.fa       # Abs path needed for Trinotate

kallisto_idx="$dir_kallisto"/"$assembly_id".idx
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
minq=5               # Min. qual for TrimGalore trimming - follows https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
minlen=36            # Min. read length for TrimGalore trimming - follows https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.htmlq
mapmax_fun=100       # Max. nr of multimapping for fungal (pre-)mapping

config_trinotate=$PWD/workflow/config/trinotate_conf.txt

# PREPROCESSING ----------------------------------------------------------------
## Run FastQC
for fastq_file in "$dir_fq"/*fastq.gz; do
    sbatch "$scr_mqc" "$fastq_file" "$dir_fqc"/raw
done

## Run MultiQC
sbatch "$scr_mqc" "$dir_fqc" "$dir_multiqc"

## Error correction with rcorrector
sbatch scripts/rcorrector.sh "$dir_fq" "$dir_rcorr"

## Remove reads with unfixable errors, according to rcorrector
for R1 in "$dir_rcorr"/*R1*.cor.fq.gz; do
    sbatch scripts/rcorrfilter.sh "$R1" "$dir_readfilt"
done

## Run TrimGalore
for R1 in "$dir_readfilt"/*R1*fastq.gz; do
    sbatch "$scr_trimg" -i "$R1" -o "$dir_trim"/rcorr -O "$dir_fqc"/trimmed_rcorr -q "$minq" -l "$minlen"
done
for R1 in "$dir_fq"/*R1*fastq.gz; do
    sbatch "$scr_trimg" -i "$R1" -o "$dir_trim"/nocorr -O "$dir_fqc"/trimmed_nocorr -q "$minq" -l "$minlen"
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


# TRANSCRIPTOME ASSEMBLY ------------------------------------------------
## Assemble with Trinity
sbatch "$bin"/trinity.sh "$dir_map2fun"/rcorr/unmapped "$dir_trinity_scratch"
cp "$dir_trinity_scratch"/Trinity.fasta "$dir_trinity"

## Redundant transcript removal with CD-HIT
sbatch "$bin"/cd-hit.sh $dir_trinity/Trinity.fasta "$assembly"


# ASSEMBLY QC -----------------------------------------------------------
## Run Trinity-stats
conda activate /users/PAS0471/jelmer/miniconda3/envs/trinotate-env
TrinityStats.pl $dir_trinity/Trinity.fasta > "$dir_trinity"/Trinity_assembly.metrics

## Run BUSCO
busco_db=insecta_odb10
sbatch "$bin"/busco.sh "$assembly" "$dir_busco" "$busco_db"
mv busco_*log "$dir_busco"; mv busco_downloads "$dir_busco"

## rnaQUAST
sbatch "$bin"/rnaquast.sh "$assembly" "$dir_rnaquast"

## Run Detonate (Note: also tried TransRate but couldn't get it to work)
sbatch "$bin"/detonate.sh "$assembly" "$dir_trim"/rcorr "$dir_detonate"

#? Try Bellerophon (again?)

## Quantify read support by mapping with Bowtie - https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
#> https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly
# bowtie2-build -–threads 4 $1 $assembly_prefix
# bowtie2 -p 10 -q --no-unal -k 20 -x $1 -1 $2 -2 $3  2>align_stats.txt| ${SINGULARITY_EXEC} samtools view -@10 -Sb -o bowtie2.bam
#? A typical Trinity transcriptome assembly will have the vast majority of all reads mapping back to the assembly, and ~70-80% of the mapped fragments found mapped as proper pairs 

## Count full-length transcripts?
#> https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts


# ASSEMBLY FILTERING ----------------------------------------------------
## - Another round of fungal sequence removal? See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6463014/, which used competitive BlastP
## - Can try TPM>1 filter using Detonate results, as in [Error, noise and bias in de novo transcriptome assemblies](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13156)
## - Can try BlastP evalue < 1e-5 as in [Error, noise and bias in de novo transcriptome assemblies](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13156)
## - Try DRAP? http://www.sigenae.org/drap/install.html

## Check for contamination with Kraken2
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/PlusPFP/
kraken_script=mcic-scripts/metagenomics/kraken-run.sh
sbatch "$kraken_script" -i "$dir_trinity"/Trinity.fasta -d "$kraken_db" -o results/kraken/trinity

fq=data/fastq/EMySc1_S23_L002_R1_001.fastq.gz
sbatch "$kraken_script" -i "$fq" -d "$kraken_db" -o results/kraken/fastq


# ASSEMBLY ANNOTATION ---------------------------------------------------
#? Try dammit? https://angus.readthedocs.io/en/2017/dammit_annotation.html#

## Trinotate
sbatch $bin/trinotate.sh "$assembly" "$config_trinotate" "$dir_trinotate"


# GENE AND TRANSCRIPT QUANTIFICATION -------------------------------------------
#? https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification

## Kallisto
sbatch mcic-scripts/rnaseq/kallisto_index.sh -i "$assembly" -o "$kallisto_idx"
for R1 in "$dir_trim"/nocorr/*R1*fastq.gz; do
    sbatch mcic-scripts/rnaseq/kallisto_quant.sh -i "$R1" -o "$dir_kallisto" -r "$kallisto_idx" 
done

# https://www.nature.com/articles/s41598-019-44499-3#Sec2
#> Subsequently, the expression abundance for each contig was estimated using one alignment-based and two alignment-free quantifiers,
#> namely (1) Bowtie2 (ver. 2.3.0)38 (–dpad 0–gbar 99999999–mp 1,1–np 1–score-min L,0,-0,1 -k 200–sensitive–no-mixed–no-discordant) followed by RSEM (ver. 1.2.31) (default parameters);
#> (2) Kallisto (ver. 0.43.0) (indexing with -k 31 and quantifying with default parameters); and (3) Salmon (ver. 0.8.2) (indexing with -k 31 and quantifying with default parameters).


# DIFFERENTIAL EXPRESSION ------------------------------------------------------
sbatch scripts/tximport.R
## Gene-to-transcript map: results/trinotate/GENE_TRANS_MAP

#? Try Ballgown? https://www.bioconductor.org/packages/devel/bioc/vignettes/ballgown/inst/doc/ballgown.html


#TODO
#> https://link.springer.com/protocol/10.1007%2F978-1-0716-1609-3_8#Sec6

