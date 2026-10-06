#!/usr/bin/env bash
# =============================================================================
# Stage III-B | MAG characterization
# =============================================================================
# Taxonomic, metabolic and evolutionary properties of the representative MAG
# set. Each step runs once for the whole set, not per sample.
# -----------------------------------------------------------------------------
set -euo pipefail

SAMPLES=(Getz1_up Getz1_down Getz2_up Getz2_down \
         Dotson1_up Dotson1_down Dotson2_up Dotson2_down)

MAG_DIR="dRep_output/dereplicated_genomes"
THREADS=120
GTDBTK_DB="/path/to/gtdbtk_db/release226"

# 1. MAG classification with GTDB  (GTDB-Tk, reference data R226)
#    Defaults apart from thread count and the genome file extension.
export GTDBTK_DATA_PATH=${GTDBTK_DB}

gtdbtk classify_wf \
    --genome_dir ${MAG_DIR} --out_dir gtdbtk_results \
    --cpus 80 --extension fa

# 2. Relative abundance  (CoverM)
#    Reads of every library are mapped against the representative MAG set,
#    listed as consecutive R1 R2 pairs after --coupled. These are the fastp
#    output before host removal.
#    --methods relative_abundance (default), --mapper bwa-mem2 (default is
#    minimap2-sr), --min-covered-fraction not set so the default 10% applies.
coverm genome \
    --threads 80 \
    --genome-fasta-directory ${MAG_DIR} --genome-fasta-extension fa \
    --coupled $(for S in "${SAMPLES[@]}"; do \
                   echo -n "${S}_clean_1.fastq.gz ${S}_clean_2.fastq.gz "; done) \
    --mapper bwa-mem2 \
    --methods relative_abundance \
    --output-file MAG_abundance_summary.csv

# 3. Metabolic annotation  (DRAM)
#    DRAM performs its own gene call internally (P2). Its gene identifiers are
#    separate from the Stage II contig-level ORF set (P1) but are shared with the
#    inStrain gene call (P3) below, so DRAM annotations and inStrain per-gene
#    statistics can be merged directly. Pathway completeness is evaluated on KEGG
#    KO identifiers. Defaults apart from thread count.
DRAM.py annotate \
    -i "${MAG_DIR}/*.fa" \
    -o DRAM_annotation \
    --threads 100

DRAM.py distill \
    -i DRAM_annotation/annotations.tsv \
    -o dram_distill_out \
    --trna_path DRAM_annotation/trnas.tsv \
    --rrna_path DRAM_annotation/rrnas.tsv

# 3b. DNA repair pathway completeness  (calculate_repair.py, this repository)
#
#     product.tsv from the distill step above holds two kinds of column:
#       - module and electron transport chain completeness, as the fraction of
#         steps with a matching orthologue;
#       - functional presence, TRUE only when every gene set DRAM lists for that
#         function is recovered (all sets satisfied, each by any one of its
#         alternative orthologues).
#
#     The four canonical DNA repair systems are KEGG *pathways* (ko03410 BER,
#     ko03420 NER, ko03430 MMR, ko03440 HR), not KEGG modules, so DRAM does not
#     evaluate them and they are absent from product.tsv. calculate_repair.py
#     scores them from the same annotations.tsv using curated step definitions
#     that carry the counterparts of both Bacteria and Archaea, and writes:
#
#       repair_completeness.tsv   per-genome fraction of steps recovered, plus a
#                                 binary call on the diagnostic steps
#       repair_steps.tsv          step-by-step detail with the KOs that matched
#       repair_diagnostics.txt    KOs and steps the annotation never recovered
#
#     Merging with product.tsv: the fraction columns use the same definition as
#     DRAM module completeness and can be displayed on one scale. The binary
#     columns do NOT — DRAM requires every listed gene set, the repair criterion
#     only the diagnostic steps of the pathway. A figure showing both must state
#     in its legend which criterion applies to which columns.
python calculate_repair.py \
    -i DRAM_annotation/annotations.tsv \
    -o repair_completeness

# 4. Microdiversity and dN/dS  (inStrain)
#    The representative MAGs are concatenated into one reference. Contig names
#    are prefixed with the MAG name so identifiers stay unique, and a
#    scaffold-to-bin (.stb) table records the mapping. Genes are re-called on
#    the concatenated reference so that coordinates match it (P3). These gene
#    identifiers are shared with the DRAM call (P2) above, but not with the
#    contig-level call of Stage II (P1).
#
#    dN/dS is computed over ALL predicted coding genes and therefore describes
#    genome-wide selection, not a repair-gene subset.
#
#    --database_mode is the inStrain preset for a reference holding many
#    genomes; it changes several defaults at once. Reads are the fastp output
#    before host removal.
mkdir -p instrain && cd instrain

: > mags.stb
: > all_mags.fasta
for F in ../${MAG_DIR}/*.fa; do
    BIN=$(basename "${F}" .fa)
    sed "s/^>/>${BIN}_/" "${F}" >> all_mags.fasta
    grep "^>" "${F}" | sed 's/>//g' | awk -v b="${BIN}" '{print b"_"$1"\t"b}' >> mags.stb
done

bowtie2-build --threads ${THREADS} all_mags.fasta all_mags_idx
prodigal -i all_mags.fasta -d all_mags_genes.fna -m -p meta -q

for S in "${SAMPLES[@]}"; do
    bowtie2 -p ${THREADS} -x all_mags_idx \
            -1 ../${S}_clean_1.fastq.gz -2 ../${S}_clean_2.fastq.gz \
      | samtools sort -@ ${THREADS} -o mapped_${S}.bam
    samtools index -@ ${THREADS} mapped_${S}.bam

    inStrain profile mapped_${S}.bam all_mags.fasta \
        -o inStrain_Profile_${S} -p ${THREADS} \
        -g all_mags_genes.fna -s mags.stb \
        --database_mode
done

# Population comparison across libraries (popANI). Profiles whose coverage
# falls below the 1x reliability threshold of inStrain are excluded; in this
# study that was Dotson1_down.
inStrain compare \
    -i $(ls -d inStrain_Profile_* | grep -v "Dotson1_down") \
    -o inStrain_Compare_Result -p ${THREADS}

inStrain genome_wide -i inStrain_Compare_Result -s mags.stb
