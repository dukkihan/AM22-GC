#!/usr/bin/env bash
# =============================================================================
# Stage I | Read preprocessing
# =============================================================================
# Contamination screening and prokaryotic read-fraction assessment before
# assembly.
#
# In  : <sample>_R1.fastq.gz, <sample>_R2.fastq.gz
# Out : <sample>_selected_1.fastq.gz, <sample>_selected_2.fastq.gz
#       (analysis-ready reads, reused by Stages II and III)
#
# Flags annotated "default" match the tool's own default and are written out
# only for the record.
# -----------------------------------------------------------------------------
set -euo pipefail

SAMPLES=(Getz1_up Getz1_down Getz2_up Getz2_down \
         Dotson1_up Dotson1_down Dotson2_up Dotson2_down)

THREADS=64
KRAKEN2_DB="/path/to/kraken2_db"
HG38_INDEX="/path/to/bowtie2/hg38"

mkdir -p qc kraken human_removed

for S in "${SAMPLES[@]}"; do

    # 1. Fastq QC  (fastp)
    #    Adapter trimming (on by default) and low-quality read filtering.
    #    -q 20 minimum base quality (default 15)
    #    -u 30 max % of bases below -q per read (default 40)
    fastp \
        -i ${S}_R1.fastq.gz -I ${S}_R2.fastq.gz \
        -o ${S}_clean_1.fastq.gz -O ${S}_clean_2.fastq.gz \
        -h qc/${S}.fastp.html -j qc/${S}.fastp.json \
        -q 20 -u 30 -w ${THREADS}

    # 2. Read-level taxonomic classification  (Kraken2)
    #    Pre-assembly screen of read composition: prokaryotic read fraction and
    #    host contamination. --confidence was not set, so the default 0.0 applies.
    kraken2 \
        --db ${KRAKEN2_DB} \
        --threads ${THREADS} \
        --paired --gzip-compressed \
        --output kraken/${S}.kraken \
        --report kraken/${S}_k2_report.txt \
        ${S}_clean_1.fastq.gz ${S}_clean_2.fastq.gz

    # 3. Human contamination removal with hg38  (Bowtie2)
    #    Read pairs that fail to align concordantly to hg38 are retained; the
    #    SAM stream is discarded.
    #    --very-sensitive (default is --sensitive)
    bowtie2 \
        -x ${HG38_INDEX} \
        -1 ${S}_clean_1.fastq.gz -2 ${S}_clean_2.fastq.gz \
        --un-conc-gz human_removed/${S}_selected_%.fastq.gz \
        --threads ${THREADS} \
        --very-sensitive \
        -S /dev/null 2> human_removed/${S}_bowtie2.log

done
