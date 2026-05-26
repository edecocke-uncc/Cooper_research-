#!/usr/bin/env python3
"""
aggregate_pan_transcriptome.py
-------------------------------
Aggregate per-accession CDS FASTAs into a single pan-transcriptome FASTA.
For each pan-gene locus, keeps the longest CDS sequence across all accessions.

Maps transcript IDs → gene IDs using the actual GFF Parent relationships
instead of regex stripping, which handles all edge cases correctly.

Usage:
    python aggregate_pan_transcriptome.py
"""

import re
import glob
import csv
from pathlib import Path

# ── PATHS ─────────────────────────────────────────────────────────────────────
CDS_DIR   = "/users/edecocke/step5/cds_per_accession"
OUT_FASTA = "/users/edecocke/step5/pan_transcriptome.fasta"
PAV_CSV   = "/users/edecocke/step5/PAV_13accessions.csv"

GFFS = [
    "/projects/cooper_research/Ware_Lab_Annotations/ChineseAmber.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/Leoti.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/pi229841.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/pi297155.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/pi300119.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/pi329311.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/pi506069.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/pi510757.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/pi655972.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ware_Lab_Annotations/Rio.pan-gene.aed_lt1.gff",
    "/projects/cooper_research/Ref_Genomes/pantranscriptome/to_Liz_sorghum_18_genomes_pangene_annotations/tx430.pan-gene.coding.gff",
    "/users/edecocke/step5/btx623_chr_only.gff",
    "/users/edecocke/step5/Grassl.pan-gene.liftoff.gff",
]
# ─────────────────────────────────────────────────────────────────────────────

ID_RE     = re.compile(r'(?:^|;)ID=([^;]+)')
PARENT_RE = re.compile(r'(?:^|;)Parent=([^;]+)')


def build_transcript_to_gene(gff_paths):
    """
    Parse GFF files and return a dict mapping transcript/mRNA ID → gene ID.
    Handles all ID formats correctly using the actual GFF Parent field.
    """
    tx2gene = {}
    for path in gff_paths:
        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                feat = parts[2]
                # Only process transcript-level features
                if feat not in ('mRNA', 'transcript'):
                    continue
                attrs = parts[8]
                id_m     = ID_RE.search(attrs)
                parent_m = PARENT_RE.search(attrs)
                if id_m and parent_m:
                    tx_id   = id_m.group(1)
                    gene_id = parent_m.group(1)
                    tx2gene[tx_id] = gene_id
    return tx2gene


def parse_fasta(path):
    """Parse a FASTA file into a dict of {seq_id: sequence}."""
    seqs = {}
    header = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if header is not None:
        seqs[header] = ''.join(buf)
    return seqs


# ── BUILD TRANSCRIPT → GENE MAPPING ──────────────────────────────────────────
print("Building transcript→gene map from GFFs ...")
tx2gene = build_transcript_to_gene(GFFS)
print(f"  {len(tx2gene):,} transcript→gene mappings loaded\n")

# ── PROCESS CDS FASTAs ────────────────────────────────────────────────────────
best        = {}   # gene_id -> longest sequence
unmapped    = {}   # tx_id -> seq  (no GFF mapping found, kept separately)

fasta_files = sorted(glob.glob(f"{CDS_DIR}/*.cds.fasta"))
if not fasta_files:
    raise FileNotFoundError(f"No *.cds.fasta files found in {CDS_DIR}")

print(f"Found {len(fasta_files)} CDS FASTA files:\n")

for fasta_path in fasta_files:
    accession = Path(fasta_path).stem.replace('.cds', '')
    seqs = parse_fasta(fasta_path)
    print(f"  {accession:15s}: {len(seqs):>7,} sequences")

    for tx_id, seq in seqs.items():
        gene_id = tx2gene.get(tx_id)
        if gene_id is None:
            unmapped[tx_id] = seq
            continue
        if gene_id not in best or len(seq) > len(best[gene_id]):
            best[gene_id] = seq

print(f"\n  Unmapped transcript IDs (not found in any GFF): {len(unmapped):,}")
if unmapped:
    print("  First 10 unmapped:")
    for tx in list(unmapped.keys())[:10]:
        print(f"    {tx}")

print(f"\nTotal unique gene loci with CDS: {len(best):,}")

# ── WRITE OUTPUT ──────────────────────────────────────────────────────────────
print(f"\nWriting to {OUT_FASTA} ...")

with open(OUT_FASTA, 'w') as out:
    for gene_id, seq in sorted(best.items()):
        out.write(f">{gene_id}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")

print("Done.")

# ── VERIFY AGAINST PAV MATRIX ─────────────────────────────────────────────────
print(f"\nVerifying against PAV matrix: {PAV_CSV}")

try:
    pav_loci = set()
    with open(PAV_CSV) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row:
                pav_loci.add(row[0])

    fasta_loci = set(best.keys())
    missing = pav_loci - fasta_loci
    extra   = fasta_loci - pav_loci

    print(f"  PAV loci:              {len(pav_loci):>7,}")
    print(f"  FASTA loci:            {len(fasta_loci):>7,}")
    print(f"  Missing from FASTA:    {len(missing):>7,}  (in PAV but no CDS extracted)")
    print(f"  Extra in FASTA:        {len(extra):>7,}  (gene ID not in PAV)")

    if missing:
        print(f"\n  First 10 missing:")
        for locus in sorted(missing)[:10]:
            print(f"    {locus}")
    if extra:
        print(f"\n  First 10 extra:")
        for locus in sorted(extra)[:10]:
            print(f"    {locus}")

    if len(missing) == 0 and len(extra) == 0:
        print("\n  ✅ FASTA loci match PAV matrix exactly.")
    else:
        print("\n  ⚠️  Mismatch — review above.")

except FileNotFoundError:
    print(f"  ⚠️  PAV file not found at {PAV_CSV} — skipping verification.")
