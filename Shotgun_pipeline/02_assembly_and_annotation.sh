#!/usr/bin/env bash
# =============================================================================
# Stage II | Assembly and contig-level annotation
# =============================================================================
# A single gene-prediction run supplies the ORF set shared by all three
# functional annotations, so KEGG, CAZy and N-cycle results merge at ORF level.
# -----------------------------------------------------------------------------
set -euo pipefail

SAMPLES=(Getz1_up Getz1_down Getz2_up Getz2_down \
         Dotson1_up Dotson1_down Dotson2_up Dotson2_down)

THREADS=80
MEMORY=680000000000                        # bytes available to MEGAHIT
KOFAM_PROFILES="/path/to/kofam/profiles"
KOFAM_KO_LIST="/path/to/kofam/ko_list"
DBCAN_HMM="/path/to/dbCAN-HMMdb-V14.txt"
NCYC_DB="/path/to/ncyc_100_db.dmnd"

mkdir -p megahit_out prodigal kofam_results dbcan_results ncyc_results

# 1. Metagenome assembly  (MEGAHIT)
#    Individual assembly per sample; no co-assembly.
#    --k-list matches MEGAHIT's default and is written out for the record.
#    --min-contig-len was not set, so the default of 200 bp applies.
for S in "${SAMPLES[@]}"; do
    megahit \
        -1 human_removed/${S}_selected_1.fastq.gz \
        -2 human_removed/${S}_selected_2.fastq.gz \
        -o megahit_out/${S}_assembly --out-prefix ${S} \
        -m ${MEMORY} -t 100 \
        --k-list 21,29,39,59,79,99,119,141        # default
done

# 2. Assembly quality assessment  (QUAST)
#    Run once over all eight assemblies for a cross-sample comparison. All
#    thresholds at default (minimum contig length 500 bp). Terminal step.
quast -o quast_results --threads 40 \
    $(for S in "${SAMPLES[@]}"; do echo -n "megahit_out/${S}_assembly/${S}.contigs.fa "; done)

for S in "${SAMPLES[@]}"; do

    # 3. Gene prediction  (Prodigal)
    #    First of the three gene calls in this study (P1); contig coordinate
    #    space. The protein set is shared as-is by steps 4-6, but its gene
    #    identifiers are not shared with the MAG-level calls in script 04.
    #    -p meta (default is -p single)
    prodigal \
        -i megahit_out/${S}_assembly/${S}.contigs.fa \
        -a prodigal/${S}_proteins.faa \
        -d prodigal/${S}_genes.fna \
        -o prodigal/${S}_prodigal.out \
        -p meta

    # 4. KEGG annotation  (KofamScan)
    #    The adaptive per-KO score thresholds stay active; -E adds a fixed
    #    E-value ceiling on top of them.
    #    --format detail-tsv (default is "detail")
    exec_annotation \
        -o kofam_results/${S}_kofam.tsv \
        -p ${KOFAM_PROFILES} -k ${KOFAM_KO_LIST} \
        --cpu ${THREADS} \
        --format detail-tsv \
        -E 1e-10 \
        prodigal/${S}_proteins.faa

    # 5. CAZy annotation  (dbCAN2 HMM database, searched with HMMER)
    #    -E 1e-5 (default is 10.0)
    hmmsearch \
        --cpu ${THREADS} \
        --domtblout dbcan_results/${S}_dbCAN.domtbl \
        --noali -E 1e-5 \
        ${DBCAN_HMM} prodigal/${S}_proteins.faa > /dev/null

    # 6. Nitrogen-cycle annotation  (NCycDB, searched with DIAMOND)
    #    --evalue 1e-10 (default 0.001), --max-target-seqs 1 (default 25),
    #    --sensitive (default is fast mode), --outfmt 6 (default)
    diamond blastp \
        --threads ${THREADS} \
        --db ${NCYC_DB} \
        --query prodigal/${S}_proteins.faa \
        --out ncyc_results/${S}_ncyc.blast6 \
        --evalue 1e-10 --max-target-seqs 1 --sensitive --outfmt 6

done
