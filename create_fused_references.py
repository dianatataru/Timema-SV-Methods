#!/usr/bin/env python3
"""
Create three versions of the HGS2 reference genome with Chr3 and Chr4 fused,
matching the fusion arrangements in RGS2, RGUS1, and RGUS2.

Fusion arrangements:
  RGUS1 version: Chr3 + revcomp(Chr4)  (Chr3 first, Chr4 reverse complemented appended)
  RGUS2 version: Chr4 + revcomp(Chr3)  (Chr4 first, Chr3 reverse complemented appended)
  RGS2  version: Chr4 + revcomp(Chr3)  (Chr4 first, Chr3 reverse complemented appended)

The fused scaffold is named Chr3_4.
"""

import os

INPUT_FASTA = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_hwy154_cen4119_hap2.fasta.masked"
OUTPUT_DIR  = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed"

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]

def parse_fasta(filepath):
    """Read fasta into ordered list of (header, sequence) tuples."""
    records = []
    header = None
    seq_parts = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
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
records = parse_fasta(INPUT_FASTA)

# Extract Chr3 and Chr4, keep everything else in order
chr3_seq = None
chr4_seq = None
other_records = []

for header, seq in records:
    name = header.split()[0]
    if name == "Chr3":
        chr3_seq = seq
        print(f"  Found Chr3: {len(seq):,} bp")
    elif name == "Chr4":
        chr4_seq = seq
        print(f"  Found Chr4: {len(seq):,} bp")
    else:
        other_records.append((header, seq))

if chr3_seq is None or chr4_seq is None:
    raise ValueError("Could not find Chr3 and/or Chr4 in input fasta. Check header names.")

# Define the three fusion versions
versions = {
    "RGUS1": (chr3_seq + revcomp(chr4_seq),         "Chr3 + revcomp(Chr4)"),
    "RGUS2": (chr4_seq + revcomp(chr3_seq),         "Chr4 + revcomp(Chr3)"),
    "RGS2":  (chr4_seq + revcomp(chr3_seq),"Chr4 + revcomp(Chr3)"),
}

for genome_name, (fused_seq, description) in versions.items():
    # Build record list: replace Chr3 with fused sequence, drop Chr4
    output_records = []
    for header, seq in records:
        name = header.split()[0]
        if name == "Chr3":
            output_records.append((f"Chr3_4", fused_seq))
            print(f"\n  {genome_name}: fusing {description} -> {len(fused_seq):,} bp")
        elif name == "Chr4":
            pass  # absorbed into Chr3
        else:
            output_records.append((header, seq))

    outfile = os.path.join(OUTPUT_DIR, f"t_crist_hwy154_cen4119_hap2.fasta.masked.fused_{genome_name}")
    write_fasta(output_records, outfile)

print("\nDone. Created 3 fused reference versions:")
