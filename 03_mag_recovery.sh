#!/usr/bin/env bash
# =============================================================================
# Stage III-A | MAG recovery
# =============================================================================
# Population genomes are recovered from contigs, then reduced to a
# non-redundant representative set.
#
# The binning signal is the abundance contrast between the two horizons of the
# same core, so reads of both horizons are mapped onto the contigs of each
# sample. CORES defines the pairing.
# -----------------------------------------------------------------------------
set -euo pipefail

CORES=(Getz1 Getz2 Dotson1 Dotson2)
THREADS=120

mkdir -p mapping binning clean_bins bin_refinement

for CORE in "${CORES[@]}"; do
  for H in up down; do

    S="${CORE}_${H}"
    [ "${H}" == "up" ] && P="${CORE}_down" || P="${CORE}_up"
    CONTIGS="megahit_out/${S}_assembly/${S}.contigs.fa"
    OUT="binning/${S}"
    mkdir -p "${OUT}/metabat2_bins" "${OUT}/maxbin2_bins"

    # 1. Read crossmapping  (Bowtie2)
    #    Reads of both horizons are mapped onto the contigs of this sample.
    #    --sensitive is the end-to-end default, written out for the record.
    bowtie2-build --threads ${THREADS} "${CONTIGS}" mapping/${S}_index

    for R in "${S}" "${P}"; do
        bowtie2 -x mapping/${S}_index \
            -1 human_removed/${R}_selected_1.fastq.gz \
            -2 human_removed/${R}_selected_2.fastq.gz \
            --threads ${THREADS} --sensitive \
          | samtools view -@ 20 -bS - \
          | samtools sort -@ 20 -o mapping/${R}_to_${S}.sorted.bam
        samtools index mapping/${R}_to_${S}.sorted.bam
    done

    # 2. Contig binning  (MetaBAT2 and MaxBin2)
    #    Two algorithms on the same coverage profile. Both at defaults apart
    #    from thread count, so MetaBAT2's 2500 bp minimum contig length applies.
    jgi_summarize_bam_contig_depths --outputDepth ${OUT}/depth.txt \
        mapping/${S}_to_${S}.sorted.bam mapping/${P}_to_${S}.sorted.bam

    metabat2 -i "${CONTIGS}" -a ${OUT}/depth.txt \
        -o ${OUT}/metabat2_bins/bin -t ${THREADS}

    # MaxBin2 takes contig ID plus the mean-depth column of each library
    cut -f1,4,6 ${OUT}/depth.txt > ${OUT}/maxbin_abund.txt
    run_MaxBin.pl -contig "${CONTIGS}" -abund ${OUT}/maxbin_abund.txt \
        -out ${OUT}/maxbin2_bins/bin -thread ${THREADS}

    # 3. Bin refinement and quality assessment  (metaWRAP, CheckM)
    #    The two bin sets are consolidated and screened on completeness and
    #    contamination. Bins passing are treated as MAGs from here on.
    #    -c 50 minimum completeness (default 70; MIMAG medium-quality draft)
    #    -x 10 maximum contamination (default)
    #    metaWRAP derives sample names from directory paths, so bin directories
    #    are named without underscores.
    TAG="${S//_/}"
    mkdir -p clean_bins/${TAG}/A clean_bins/${TAG}/B
    cp ${OUT}/metabat2_bins/*.fa   clean_bins/${TAG}/A/ 2>/dev/null || true
    cp ${OUT}/maxbin2_bins/*.fasta clean_bins/${TAG}/B/ 2>/dev/null || true

    metawrap bin_refinement \
        -o bin_refinement/${TAG} -t ${THREADS} -m 700 \
        -A clean_bins/${TAG}/A/ -B clean_bins/${TAG}/B/ \
        -c 50 -x 10

  done
done

# 4. Representative MAG selection with duplicate removal  (dRep)
#    Run once on the MAGs of all eight samples pooled. File names are prefixed
#    with the sample so that bin names stay unique.
#    -sa 0.95 secondary clustering ANI, species level (default)
#    -nc 0.30 minimum overlap between genomes (default 0.10)
#
#    Completeness and contamination filters were left at their defaults
#    (-comp 75, -con 25). Refinement admits MAGs at >= 50% completeness, so it
#    is dRep's default completeness filter that sets the >= 75% completeness of
#    the final representative set.
#
#    genome_quality.csv is the CheckM table from the refinement step, in dRep
#    genomeInfo format: genome,completeness,contamination
mkdir -p dRep_input
for CORE in "${CORES[@]}"; do
  for H in up down; do
    TAG="${CORE}${H}"
    for B in bin_refinement/${TAG}/metawrap_50_10_bins/*.fa; do
        cp "${B}" dRep_input/${TAG}_$(basename "${B}")
    done
  done
done

dRep dereplicate dRep_output \
    -g dRep_input/*.fa \
    --genomeInfo genome_quality.csv \
    -sa 0.95 -nc 0.30 -p ${THREADS}
