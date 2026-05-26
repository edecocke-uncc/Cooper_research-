#!/usr/bin/env python3
"""
build_pav_from_gffs.py
-----------------------
Builds a presence/absence variation (PAV) matrix directly from per-accession
GFF files. For each GFF, extracts all gene-level IDs. Presence = gene ID found
in that accession's GFF; absence = not found.

Output: /users/edecocke/step5/PAV_13accessions.csv

Rows    = union of all gene locus IDs across all 13 GFFs
Columns = one per accession (0/1)
"""

import re
import csv
from collections import defaultdict

# ── GFF PATHS ─────────────────────────────────────────────────────────────────
GFFS = {
    "ChineseAmber": "/projects/cooper_research/Ware_Lab_Annotations/ChineseAmber.pan-gene.aed_lt1.gff",
    "Leoti":        "/projects/cooper_research/Ware_Lab_Annotations/Leoti.pan-gene.aed_lt1.gff",
    "pi229841":     "/projects/cooper_research/Ware_Lab_Annotations/pi229841.pan-gene.aed_lt1.gff",
    "pi297155":     "/projects/cooper_research/Ware_Lab_Annotations/pi297155.pan-gene.aed_lt1.gff",
    "pi300119":     "/projects/cooper_research/Ware_Lab_Annotations/pi300119.pan-gene.aed_lt1.gff",
    "pi329311":     "/projects/cooper_research/Ware_Lab_Annotations/pi329311.pan-gene.aed_lt1.gff",
    "pi506069":     "/projects/cooper_research/Ware_Lab_Annotations/pi506069.pan-gene.aed_lt1.gff",
    "pi510757":     "/projects/cooper_research/Ware_Lab_Annotations/pi510757.pan-gene.aed_lt1.gff",
    "pi655972":     "/projects/cooper_research/Ware_Lab_Annotations/pi655972.pan-gene.aed_lt1.gff",
    "Rio":          "/projects/cooper_research/Ware_Lab_Annotations/Rio.pan-gene.aed_lt1.gff",
    "tx430":        "/projects/cooper_research/Ref_Genomes/pantranscriptome/to_Liz_sorghum_18_genomes_pangene_annotations/tx430.pan-gene.coding.gff",
    "btx623":       "/users/edecocke/step5/btx623_chr_only.gff",
    "Grassl":       "/users/edecocke/step5/Grassl.pan-gene.liftoff.gff",
}

OUT_CSV = "/users/edecocke/step5/PAV_13accessions.csv"

# Regex to extract ID= from GFF attributes column
ID_RE = re.compile(r'(?:^|;)ID=([^;]+)')


def get_gene_ids(gff_path):
    """Return set of gene IDs from a GFF (lines where feature type == 'gene')."""
    gene_ids = set()
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'gene':
                continue
            m = ID_RE.search(parts[8])
            if m:
                gene_ids.add(m.group(1))
    return gene_ids


# ── MAIN ──────────────────────────────────────────────────────────────────────

accessions = list(GFFS.keys())
presence = {}   # accession -> set of gene IDs

print("Parsing GFFs:\n")
for acc, path in GFFS.items():
    ids = get_gene_ids(path)
    presence[acc] = ids
    print(f"  {acc:15s}: {len(ids):>7,} gene loci")

# Union of all gene IDs
all_loci = set()
for ids in presence.values():
    all_loci.update(ids)

print(f"\nTotal unique loci across all 13 accessions: {len(all_loci):,}")

# ── WRITE PAV CSV ─────────────────────────────────────────────────────────────
print(f"\nWriting PAV to {OUT_CSV} ...")

with open(OUT_CSV, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Gene'] + accessions)
    for locus in sorted(all_loci):
        row = [locus] + [1 if locus in presence[acc] else 0 for acc in accessions]
        writer.writerow(row)

print("Done.")

# ── SUMMARY STATS ─────────────────────────────────────────────────────────────
import pandas as pd
df = pd.read_csv(OUT_CSV)
acc_cols = df.columns[1:]
df['n_present'] = df[acc_cols].sum(axis=1)

core       = (df['n_present'] == len(accessions)).sum()
accessory  = ((df['n_present'] > 0) & (df['n_present'] < len(accessions))).sum()
private    = (df['n_present'] == 1).sum()

print(f"\nPAV summary:")
print(f"  Total loci:       {len(df):>7,}")
print(f"  Core (all 13):    {core:>7,}")
print(f"  Accessory:        {accessory:>7,}")
print(f"  Private (1 acc):  {private:>7,}")
print(f"\nPer-accession gene counts:")
for acc in accessions:
    n = df[acc].sum()
    print(f"  {acc:15s}: {n:>7,}")
