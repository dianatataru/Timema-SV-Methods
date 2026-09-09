#!/usr/bin/env python

import argparse
import pandas as pd
from collections import defaultdict
from typing import Union
import gzip

def read_vcf_line_by_line(vcf_path):
    open_func = gzip.open if vcf_path.endswith(".gz") else open
    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            yield {
                "CHROM": fields[0],
                "POS": int(fields[1]),
                "ID": fields[2],
                "REF": fields[3],
                "ALT": fields[4],
                "QUAL": fields[5],
                "FILTER": fields[6],
                "INFO": fields[7],
            }

# Part 1: Summary for all variant edges
def variant_edges_summary_from_dict(var_list: list, var_dict: dict):
    summary_dict = dict()
    for edge in sorted(list(var_list)):
        summary_dict[var_dict[edge]] = summary_dict.get(var_dict[edge], 0) + 1
    summary_dict['Total'] = len(var_list)
    return summary_dict

def prepare_dataframe_dict(var_dict):
    return pd.DataFrame({
        "Variant Type": ['SNP', 'MNP', 'Insertion', 'Deletion', 'Replacement', 'Inversion', 'Duplication', 'Total'],
        "Count": [
                  var_dict.get('SNP', 0),
                  var_dict.get('MNP', 0),
                  var_dict.get('INS', 0),
                  var_dict.get('DEL', 0),
                  var_dict.get('REP', 0),
                  var_dict.get('INV', 0),
                  var_dict.get('DUP', 0),
                  var_dict.get('Total', 0),
        ]
    })

def comprehensive_summary(graph_vcf_input: Union[str, pd.DataFrame], 
                          chrom_label: str, 
                          tandem_repeat: bool=False):
    variant_sets = defaultdict(set)
    var_dict = dict()

    if isinstance(graph_vcf_input, str):
        row_iter = read_vcf_line_by_line(graph_vcf_input)
    elif isinstance(graph_vcf_input, pd.DataFrame):
        row_iter = graph_vcf_input.to_dict(orient="records")
    else:
        raise ValueError("Input must be either a path to VCF or a pandas DataFrame.")

    for row in row_iter:
        edge = row['ID']
        info_dict = {k: v for k, v in (field.split('=') for field in row['INFO'].split(';') if '=' in field)}
        var_type = info_dict['VT']
        nearly_identical = int(info_dict['NIA'])
        allele_count = int(info_dict['AC']) if var_type == 'INV' else min(int(info_dict['RC']), int(info_dict['AC']))
        ref_allele = row['REF'] if row['REF'] != '.' else info_dict['NR']
        alt_allele = row['ALT']
        allele_length = len(ref_allele) + len(alt_allele)

        REF_PREFIX = "Hap2_t_crist_hwy154_cen4119"
        ref_cols = [col for col in row.keys() if col.startswith(REF_PREFIX + "#")]

        on_linear = False
        for col in ref_cols:
            gt = row[col].split(":")[0]   # FORMAT is GT:CR:CA
            if gt == "1":
                on_linear = True
                break

        is_repeat = info_dict['TR_MOTIF'] != '.'

        if tandem_repeat and not is_repeat:
            continue

        var_dict[edge] = var_type
        variant_sets['All'].add(edge)
        variant_sets['Linear' if on_linear else 'Off_Linear'].add(edge)
        variant_sets['Small' if allele_length < 50 or nearly_identical else 'Large'].add(edge)
        variant_sets['Common' if allele_count >= 5 else 'Uncommon'].add(edge)

    # Generate combinations
    for a in ['Linear', 'Off_Linear']:
        for b in ['Small', 'Large']:
            variant_sets[f'{a}_{b}'] = variant_sets[a].intersection(variant_sets[b])
        for c in ['Common', 'Uncommon']:
            variant_sets[f'{a}_{c}'] = variant_sets[a].intersection(variant_sets[c])
    for b in ['Small', 'Large']:
        for c in ['Common', 'Uncommon']:
            variant_sets[f'{b}_{c}'] = variant_sets[b].intersection(variant_sets[c])
            for a in ['Linear', 'Off_Linear']:
                variant_sets[f'{a}_{b}_{c}'] = variant_sets[a].intersection(variant_sets[f'{b}_{c}'])

    # Construct DataFrames
    all_var_df = prepare_dataframe_dict(variant_edges_summary_from_dict(variant_sets['All'], var_dict))
    summary_dict = {'CHROM': [chrom_label] * len(all_var_df), 'Variant Type': all_var_df['Variant Type']}
    for key, variant_set in variant_sets.items():
        df = prepare_dataframe_dict(variant_edges_summary_from_dict(variant_set, var_dict))
        summary_dict[key+"_Variants"] = df['Count']

    return pd.DataFrame(summary_dict, columns=[
        "CHROM",
        "Variant Type",
        "All_Variants",
        "Linear_Variants",
        "Off_Linear_Variants",
        "Small_Variants",
        "Large_Variants",
        "Common_Variants",
        "Uncommon_Variants",
        "Linear_Small_Variants",
        "Linear_Large_Variants",
        "Off_Linear_Small_Variants",
        "Off_Linear_Large_Variants",
        "Linear_Common_Variants",
        "Linear_Uncommon_Variants",
        "Off_Linear_Common_Variants",
        "Off_Linear_Uncommon_Variants",
        "Small_Common_Variants",
        "Small_Uncommon_Variants",
        "Large_Common_Variants",
        "Large_Uncommon_Variants",
        "Linear_Small_Common_Variants",
        "Linear_Large_Common_Variants",
        "Linear_Small_Uncommon_Variants",
        "Linear_Large_Uncommon_Variants",
        "Off_Linear_Small_Common_Variants",
        "Off_Linear_Large_Common_Variants",
        "Off_Linear_Small_Uncommon_Variants",
        "Off_Linear_Large_Uncommon_Variants"
        ])

def main():
    parser = argparse.ArgumentParser(description="Summarize Pantree SV VCF")
    parser.add_argument("--vcf", required=True, help="Path to pantree VCF or VCF.GZ")
    parser.add_argument("--chrom", required=True, help="Chromosome label (e.g. chr5)")
    parser.add_argument("--outdir", default="summary", help="Output directory")

    args = parser.parse_args()

    vcf_path = args.vcf
    chrom_label = args.chrom
    outdir = args.outdir

    import os
    os.makedirs(outdir, exist_ok=True)

    print(f"Reading VCF: {vcf_path}")
    print(f"Chrom label: {chrom_label}")
    print(f"Output dir: {outdir}")

    df = comprehensive_summary(vcf_path, chrom_label)

    outfile = f"{outdir}/pantree_sv_summary_{chrom_label}.tsv"
    df.to_csv(outfile, sep="\t", index=False)

    print(f"Written: {outfile}")

if __name__ == "__main__":
    main()
