#!/usr/bin/env python3
FASTA_IN = "/users/edecocke/step5/pan_transcriptome.fasta"
GTF_OUT  = "/users/edecocke/step5/pseudo-annotation.gtf"

seqs = {}
header = None
buf = []

with open(FASTA_IN) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if header is not None:
                seqs[header] = len(''.join(buf))
            header = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if header is not None:
        seqs[header] = len(''.join(buf))

with open(GTF_OUT, 'w') as out:
    for i, (gene_id, length) in enumerate(seqs.items(), start=1):
        attr = f'transcript_id "{gene_id}"; gene_id "{gene_id}"; gene_name "{gene_id}";'
        out.write(f"pangene{i}\tpseudo\tgene\t1\t{length}\t.\t+\t.\t{attr}\n")
        out.write(f"pangene{i}\tpseudo\texon\t1\t{length}\t.\t+\t.\t{attr}\n")

print(f"Done — {len(seqs):,} gene+exon pairs written to {GTF_OUT}")
