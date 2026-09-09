#!/usr/bin/env python3
"""
Reverse complement specific chromosomes in a fasta file.
Used to fix strand orientation issues flagged by SyRI.
"""

import os

INPUT_FASTA = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_refug_cen4120_hap1.fasta.masked"
OUTPUT_FASTA = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_refug_cen4120_hap1.fasta.masked.reoriented"

# Chromosomes to reverse complement (as they appear in the renamed fasta headers)
TO_REVCOMP = {"Chr6", "Chr8", "Chr3_4"}

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]

def parse_fasta(filepath):
    """Read fasta into list of (header, sequence) tuples, preserving order."""
    records = []
    header = None
    seq_parts = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]  # strip ">"
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            records.append((header, "".join(seq_parts)))
    return records

def write_fasta(records, filepath, line_width=50):
    with open(filepath, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + "\n")

print("Reading input fasta...")
records = parse_fasta(INPUT_FASTA)
print(f"  Loaded {len(records)} chromosomes")

updated = []
for header, seq in records:
    chr_name = header.split()[0]  # handle any trailing info in header
    if chr_name in TO_REVCOMP:
        print(f"  Reverse complementing {chr_name} ({len(seq):,} bp)")
        seq = revcomp(seq)
    else:
        print(f"  Keeping {chr_name} as-is ({len(seq):,} bp)")
    updated.append((header, seq))

print(f"\nWriting output to {OUTPUT_FASTA}...")
write_fasta(updated, OUTPUT_FASTA)
print("Done.")
