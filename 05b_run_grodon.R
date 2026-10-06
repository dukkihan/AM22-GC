# =============================================================================
# Stage III-C, step 2 | Maximum growth rate prediction  (gRodon)
# =============================================================================
# gRodon infers the maximum growth rate a genome is capable of from the codon
# usage bias of its ribosomal protein genes, relative to the rest of its coding
# sequences. Inputs are produced by 05a_extract_cds_and_ribosomal_genes.py.
#
# Reference genomes retrieved from NCBI are processed identically to the MAGs:
# cultured representatives of the recovered lineages provide an external reference,
# and Escherichia coli K-12 serves as the normalising genome.
#
# This script outputs predicted doubling times. The growth rate ratios reported in
# the paper are these values normalised to Escherichia coli K-12.
#
# predictGrowth mode:
#   "full"    complete reference genomes            (gRodon's default)
#   "partial" incomplete genomes; used for the MAGs, which are drafts
#
# Genomes with fewer than five ribosomal protein genes are skipped, since the
# codon usage estimate is unreliable below that.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(gRodon)
    library(Biostrings)
})

work_dir  <- "/path/to/gRodon"
input_dir <- file.path(work_dir, "extracted_inputs")
out_csv   <- file.path(work_dir, "gRodon_results.csv")

cds_files <- list.files(input_dir, pattern = "_cds\\.fasta$", full.names = TRUE)
results_list <- list()

for (cds_path in sort(cds_files)) {
    genome_id <- gsub("_cds\\.fasta$", "", basename(cds_path))
    ribo_path <- file.path(input_dir, paste0(genome_id, "_ribo_ids.txt"))

    if (!file.exists(ribo_path)) next

    genes    <- readDNAStringSet(cds_path)
    ribo_ids <- readLines(ribo_path)
    ribo_ids <- ribo_ids[nchar(trimws(ribo_ids)) > 0]

    gene_names <- sapply(strsplit(names(genes), "\\s+"), `[`, 1)

    # highly_expressed must be a logical vector over the gene set
    is_highly_expressed <- gene_names %in% ribo_ids
    ribo_count <- sum(is_highly_expressed)

    if (ribo_count < 5) {
        cat(sprintf("[SKIP] %-32s : Too few ribosomal genes (%d)\n", genome_id, ribo_count))
        next
    }

    # complete reference genomes are run in full mode, draft MAGs in partial mode
    mode_type <- if (grepl("^Ref_", genome_id)) "full" else "partial"

    tryCatch({
        pred <- predictGrowth(genes = genes,
                              highly_expressed = is_highly_expressed,
                              mode = mode_type)

        results_list[[genome_id]] <- data.frame(
            Genome           = genome_id,
            Mode             = mode_type,
            Doubling_Time_hr = round(pred$d, 2),
            Lower_CI_hr      = round(pred$LowerCI, 2),
            Upper_CI_hr      = round(pred$UpperCI, 2),
            CUB              = round(pred$CUB, 4),
            Ribosomal_Count  = ribo_count,
            Total_CDS_Count  = length(genes),
            stringsAsFactors = FALSE
        )

        cat(sprintf("[OK] %-32s | %-7s | d_max: %6.2f hr | Ribo: %2d\n",
                    genome_id, mode_type, pred$d, ribo_count))
    }, error = function(e) {
        cat(sprintf("[FAIL] %-32s : %s\n", genome_id, e$message))
    })
}

if (length(results_list) > 0) {
    final_df <- do.call(rbind, results_list)
    final_df <- final_df[order(final_df$Doubling_Time_hr), ]
    write.csv(final_df, out_csv, row.names = FALSE)
}
