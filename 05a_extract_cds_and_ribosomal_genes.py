#!/usr/bin/env python3
# =============================================================================
# Stage III-C, step 1 | CDS and ribosomal-protein gene extraction
# =============================================================================
# Prepares the two inputs gRodon needs for every genome:
#   <genome>_cds.fasta     all coding sequences
#   <genome>_ribo_ids.txt  identifiers of the ribosomal protein genes, used as
#                          the highly expressed gene set
#
# Two sources are handled so that the recovered MAGs and the external control
# genomes are processed identically:
#   - MAGs : coding sequences from the DRAM gene call (genes.fna), with
#            ribosomal proteins identified from the DRAM annotation table
#   - NCBI reference genomes : CDS features extracted from GenBank files
#
# Requires biopython.
# -----------------------------------------------------------------------------
import os
import glob
import re
from Bio import SeqIO

BASE_DIR   = "/path/to/gRodon"
OUT_DIR    = os.path.join(BASE_DIR, "extracted_inputs")
DRAM_ANNOT = "/path/to/DRAM_annotation/annotations.tsv"
DRAM_FNA   = "/path/to/DRAM_annotation/genes.fna"
NCBI_DIR   = os.path.join(BASE_DIR, "ncbi_ref")

os.makedirs(OUT_DIR, exist_ok=True)

broad_ribo_pattern = re.compile(
    r'(ribosomal protein [ls]\d+|50s ribosomal protein|30s ribosomal protein|'
    r'ribosome.*protein|ribosomal subunit|polypeptide.*ribosom|ribosomal protein [a-z0-9]+)',
    re.IGNORECASE
)

mag_ribo_genes = {}
gene_to_mag = {}

# -----------------------------------------------------------------------------
# MAGs: identify ribosomal protein genes from the DRAM annotation table.
# A gene is called ribosomal if its KO identifier falls in the ribosomal protein
# range, or if its KEGG or Pfam hit text matches the pattern above without also
# matching a modifying-enzyme keyword.
# -----------------------------------------------------------------------------
if os.path.exists(DRAM_ANNOT):
    with open(DRAM_ANNOT, "r") as f:
        # rstrip("\n") rather than strip() so that the leading index column,
        # which is separated by a tab, is preserved
        header_line = f.readline().rstrip("\n")
        headers = header_line.split("\t")

        ko_idx    = headers.index("ko_id")     if "ko_id"     in headers else -1
        fasta_idx = headers.index("fasta")     if "fasta"     in headers else -1
        kegg_idx  = headers.index("kegg_hit")  if "kegg_hit"  in headers else -1
        pfam_idx  = headers.index("pfam_hits") if "pfam_hits" in headers else -1

        for line in f:
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 2:
                continue
            gene_id = toks[0]

            mag_name = toks[fasta_idx] if fasta_idx != -1 and len(toks) > fasta_idx \
                       else gene_id.split("_k141")[0]
            mag_name = mag_name.replace(".fa", "").replace(".fasta", "").replace(".bin", "")
            gene_to_mag[gene_id] = mag_name

            if mag_name not in mag_ribo_genes:
                mag_ribo_genes[mag_name] = set()

            ko_id     = toks[ko_idx]   if ko_idx   != -1 and len(toks) > ko_idx   else ""
            kegg_hit  = toks[kegg_idx] if kegg_idx != -1 and len(toks) > kegg_idx else ""
            pfam_hit  = toks[pfam_idx] if pfam_idx != -1 and len(toks) > pfam_idx else ""

            is_ribo = False
            if re.match(r'K028[6-9][0-9]|K029[0-9]{2}|K030[0-9]{2}|K01977|K01988', ko_id):
                is_ribo = True
            elif broad_ribo_pattern.search(kegg_hit) or broad_ribo_pattern.search(pfam_hit):
                if not any(x in (kegg_hit + pfam_hit).lower()
                           for x in ["methyl", "synthetase", "kinase", "assembly", "modification"]):
                    is_ribo = True

            if is_ribo:
                mag_ribo_genes[mag_name].add(gene_id)

    print(f"[INFO] Parsed DRAM annotations.tsv. Found ribosomal genes for {len(mag_ribo_genes)} MAGs.")

# -----------------------------------------------------------------------------
# MAGs: write the coding sequences. genes.fna is read directly so that gene
# identifiers match the annotation table exactly. Sequences shorter than 90 bp
# or not a multiple of three are dropped.
# -----------------------------------------------------------------------------
mag_cds_records = {}
if os.path.exists(DRAM_FNA):
    print("[INFO] Parsing DRAM genes.fna directly to ensure perfect ID matching...")
    for record in SeqIO.parse(DRAM_FNA, "fasta"):
        gene_id = record.id
        mag_name = gene_to_mag.get(gene_id, gene_id.split("_k141")[0])
        mag_name = mag_name.replace(".fa", "").replace(".fasta", "").replace(".bin", "")

        seq = str(record.seq)
        if len(seq) >= 90 and len(seq) % 3 == 0:
            mag_cds_records.setdefault(mag_name, []).append(f">{gene_id}\n{seq}")

    for mag_name, records in mag_cds_records.items():
        if len(records) > 0:
            with open(os.path.join(OUT_DIR, f"{mag_name}_cds.fasta"), "w") as f:
                f.write("\n".join(records) + "\n")
            ribos = mag_ribo_genes.get(mag_name, set())
            with open(os.path.join(OUT_DIR, f"{mag_name}_ribo_ids.txt"), "w") as f:
                f.write("\n".join(ribos) + "\n")
            print(f"[OK] {mag_name}: CDS={len(records)}, Ribosomal genes={len(ribos)}")

# -----------------------------------------------------------------------------
# NCBI reference genomes: extract CDS features from GenBank files and apply the
# same length filter and ribosomal-protein pattern to the feature qualifiers.
# -----------------------------------------------------------------------------
print("\n[INFO] Parsing NCBI Reference Genomes...")
ncbi_files = glob.glob(os.path.join(NCBI_DIR, "*.gbk"))
for gbk_path in ncbi_files:
    base_name = os.path.splitext(os.path.basename(gbk_path))[0]
    cds_records = []
    ribo_ids = set()

    for record in SeqIO.parse(gbk_path, "genbank"):
        for feat in record.features:
            if feat.type == "CDS":
                gene_id = None
                for key in ["protein_id", "locus_tag", "gene", "ID"]:
                    if key in feat.qualifiers:
                        gene_id = feat.qualifiers[key][0]
                        break
                if not gene_id:
                    continue

                try:
                    seq_dna = feat.extract(record.seq)
                    if len(seq_dna) >= 90 and len(seq_dna) % 3 == 0:
                        cds_records.append(f">{gene_id}\n{str(seq_dna)}")
                except Exception:
                    continue

                all_text = " ".join([v for k, val_list in feat.qualifiers.items() for v in val_list])
                if broad_ribo_pattern.search(all_text):
                    if not any(x in all_text.lower()
                               for x in ["methyl", "synthetase", "kinase", "assembly", "modification"]):
                        ribo_ids.add(gene_id)

    if cds_records:
        with open(os.path.join(OUT_DIR, f"{base_name}_cds.fasta"), "w") as f:
            f.write("\n".join(cds_records) + "\n")
        with open(os.path.join(OUT_DIR, f"{base_name}_ribo_ids.txt"), "w") as f:
            f.write("\n".join(ribo_ids) + "\n")
        print(f"[OK] {base_name}: CDS={len(cds_records)}, Ribosomal genes={len(ribo_ids)}")

print("\n[DONE] Extraction complete.")
