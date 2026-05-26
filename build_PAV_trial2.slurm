#!/bin/bash
#SBATCH --job-name=build_PAV
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/users/edecocke/step4/logs/build_PAV_%j.out
#SBATCH --error=/users/edecocke/step4/logs/build_PAV_%j.err

mkdir -p /users/edecocke/step4/logs

module load anaconda3/2023.09

python3 - << 'PYEOF'
import pandas as pd

MATRIX    = "/projects/cooper_research/Ref_Genomes/pantranscriptome/pangenome_annotation_matrix.txt"
BTX_GFF   = "/projects/cooper_research/Ref_Genomes/pantranscriptome/to_Liz_sorghum_18_genomes_pangene_annotations/btx623.pan-gene.coding.gff"
TX430_GFF = "/projects/cooper_research/Ref_Genomes/pantranscriptome/to_Liz_sorghum_18_genomes_pangene_annotations/tx430.pan-gene.coding.gff"
OUT       = "/users/edecocke/step4/PAV_10accessions.csv"

def parse_gene_ids(gff_path):
    ids = set()
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            for attr in parts[8].split(';'):
                if attr.startswith('ID='):
                    ids.add(attr[3:])
                    break
    return ids

print("Loading existing matrix...")
df = pd.read_csv(MATRIX, sep=r'\s+', engine='python')
print(f"  {len(df)} loci, columns: {list(df.columns)}")

print("Parsing BTx623 pan-gene GFF...")
btx_ids = parse_gene_ids(BTX_GFF)
print(f"  {len(btx_ids)} loci present in BTx623")

print("Parsing TX430 pan-gene GFF...")
tx430_ids = parse_gene_ids(TX430_GFF)
print(f"  {len(tx430_ids)} loci present in TX430")

df['btx623'] = df['Gene'].apply(lambda x: 1 if x in btx_ids else 0)
df['tx430']  = df['Gene'].apply(lambda x: 1 if x in tx430_ids else 0)

print(f"\nFinal matrix: {len(df)} loci x {len(df.columns)-1} accessions")
print(f"Columns: {list(df.columns)}")
print(f"\nPer-accession presence counts:")
for col in df.columns[1:]:
    print(f"  {col}: {df[col].sum()} / {len(df)} present")

df.to_csv(OUT, index=False)
print(f"\nDone. Written to {OUT}")
PYEOF
