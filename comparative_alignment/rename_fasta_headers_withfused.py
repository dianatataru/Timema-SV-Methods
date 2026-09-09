#!/usr/bin/env python3
"""
Rename fasta scaffold headers to homologous chromosome IDs for SyRI.
Scaffold_N__... -> >Chr1, >Chr2, etc. based on lookup table.
For genomes with Chr3/Chr4 fusion on scaffold 1, scaffold 1 is renamed to Chr3.
"""

import os
import re

# --- Configuration ---
INPUT_DIR = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes"
OUTPUT_DIR = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed"

# Chromosome mapping table: Chr -> {genome: scaffold_number}
# Columns: Chr | UGS | RGS1 | RGS2 | RGUS1 | RGUS2 | HGS1 | HGS2 | HGUS1 | HGUS2
CHROM_TABLE = [
    #Chr   UGS    RGS1  RGS2  RGUS1  RGUS2  HGS1  HGS2  HGUS1  HGUS2
    ( 1,  8483,     6,    8,   11,     9,    13,   13,    22,    15),
    ( 2, 14640,     8,    5,    7,     5,     5,    6,    23,     1),
    ( 3, 42935,     2,    1,    1,     1,     3,    2,    16,     3),
    ( 4, 42912,     1,    1,    1,     1,     1,    1,    64,    35),
    ( 5, 18722,    13,   12,   12,    12,    12,   12,     5,    10),
    ( 6,  9928,     7,    2,    6,     4,     4,    5,    11,    44),
    ( 7, 10660,    10,    7,   10,    11,    10,    8,    54,     7),
    ( 8,  7748,    11,    9,    3,     3,    11,    4,     7,    23),
    ( 9, 16151,     4,    3,    8,    10,     8,    9,    46,    21),
    (10, 14160,    12,   10,    4,     6,     7,    7,    15,    16),
    (11, 12033,     9,    6,    5,     7,     9,   10,     2,    12),
    (12, 12380,     5,    4,    9,     8,     6,   11,     1,    36),
    (13, 14101,     3,   11,    2,     2,     2,    3,    36,     8),
]

# Genome file mapping: filename -> (genome_name, column index in CHROM_TABLE)
GENOMES = {
#    "t_crist_hwy154_cen4119_hap1.fasta.masked": ("HGS1",  6),
#    "t_crist_hwy154_cen4119_hap2.fasta.masked": ("HGS2",  7),
#    "t_crist_hwy154_cen4280_hap1.fasta.masked": ("HGUS1", 8),
#    "t_crist_hwy154_cen4280_hap2.fasta.masked": ("HGUS2", 9),
    "t_crist_refug_cen4122_hap1.fasta.masked":  ("RGS1",  2),
    "t_crist_refug_cen4122_hap2.fasta.masked":  ("RGS2",  3),
    "t_crist_refug_cen4120_hap1.fasta.masked":  ("RGUS1", 4),
    "t_crist_refug_cen4120_hap2.fasta.masked":  ("RGUS2", 5),
}

# For fused genomes, scaffold 1 contains Chr3+Chr4 fused — override its name to Chr3.
# Key: genome_name, Value: dict of scaffold_number -> chr_name overrides
FUSED_OVERRIDES = {
    "RGS2":  {1: "Chr3_4"},
    "RGUS1": {1: "Chr3_4"},
    "RGUS2": {1: "Chr3_4"},
}

os.makedirs(OUTPUT_DIR, exist_ok=True)

for filename, (genome_name, col_idx) in GENOMES.items():
    input_path = os.path.join(INPUT_DIR, filename)
    output_path = os.path.join(OUTPUT_DIR, filename)

    # Build scaffold_number -> ChrN mapping for this genome from the table
    scaffold_to_chr = {}
    for row in CHROM_TABLE:
        chr_num = row[0]
        scaffold_num = row[col_idx]
        # Skip duplicate scaffold numbers (fused chrs) — handled by FUSED_OVERRIDES
        if scaffold_num not in scaffold_to_chr:
            scaffold_to_chr[scaffold_num] = f"Chr{chr_num}"

    # Apply any fusion overrides for this genome
    if genome_name in FUSED_OVERRIDES:
        for scaffold_num, chr_name in FUSED_OVERRIDES[genome_name].items():
            scaffold_to_chr[scaffold_num] = chr_name
            print(f"  [{genome_name}] Override: Scaffold_{scaffold_num} -> {chr_name} (fused Chr3+Chr4)")

    print(f"Processing {genome_name}: {filename}")

    renamed = 0
    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                match = re.match(r">Scaffold_(\d+)__", line) or \
                        re.match(r">ScIzd6S_(\d+)_HRSCAF_", line) or \
			re.match(r">ScrX45T_(\d+)_HRSCAF_", line)
                if match:
                    scaffold_num = int(match.group(1))
                    if scaffold_num in scaffold_to_chr:
                        new_header = f">{scaffold_to_chr[scaffold_num]}\n"
                        fout.write(new_header)
                        renamed += 1
                    else:
                        print(f"  WARNING: Scaffold_{scaffold_num} not found in mapping table, keeping original header")
                        fout.write(line)
                else:
                    print(f"  WARNING: Could not parse header: {line.strip()}, keeping original")
                    fout.write(line)
            else:
                fout.write(line)

    print(f"  Done: {renamed} chromosomes renamed -> {output_path}\n")

print("All genomes processed.")
print(f"Renamed fastas are in: {OUTPUT_DIR}")
