#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=summarizecoords
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizecoords%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizecoords-%j.out

module load cactus/3.0.1

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree

# Output table header
echo -e "scaffold\tvariant_id\tstart_node\tref_startpos\tstart_steps\tend_node\tref_endpos\tend_steps" > ref_inversion_coordinates.tsv

# Loop over each scaffold's inversion VCF
for vcf in *_pantree_inversions_only.vcf.gz; do

    # Extract scaffold name/number from filename
    scaffold=$(basename "$vcf" .vcf.gz | sed 's/_pantree_inversions_only//')

    # Path to corresponding xg index, vg file, and display paths file
    xg="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/${scaffold}.xg"
    vg_file="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/${scaffold}.vg"
    ref_path="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/ref_display_paths/display_paths_${scaffold}.txt"

    # Check vg file exists
    if [[ ! -f "$vg_file" ]]; then
        echo "WARNING: $vg_file not found, skipping $scaffold" >&2
        continue
    fi

    ref_path_name=$(cat "$ref_path")

    # Function to find nearest ref position with step tracking
    find_ref_pos() {
        local node=$1
        local xg=$2
        local ref_path_name=$3
        local max_steps=20

        for ((c=1; c<=max_steps; c++)); do
            pos=$(vg find -x "$xg" -n "$node" -c "$c" 2>/dev/null | \
                vg position -x "$xg" -g -p "$ref_path_name" - 2>/dev/null | \
                awk 'NR==1{print $2}')
            if [[ -n "$pos" ]]; then
                echo "${pos},${c}"  # return both position and steps needed
                return
            fi
        done
        echo "NO_REF_ANCHOR,${max_steps}+"
    }

    # Loop over each inversion in the VCF (skip header lines, decompress on the fly)
    while IFS=$'\t' read -r chrom pos id ref alt qual filter info format samples; do

        # Extract start and end nodes from the ID field e.g. >2993340<2993702
        start_node=$(echo "$id" | grep -oP '(?<=[><])\d+' | head -1)
        end_node=$(echo "$id"   | grep -oP '(?<=[><])\d+' | tail -1)

        # Try direct reference position first
	ref_startpos=$(vg find -x "$xg" -n "$start_node" -P "$ref_path_name" 2>/dev/null | awk '{print $2}')
	ref_endpos=$(vg find -x "$xg" -n "$end_node" -P "$ref_path_name" 2>/dev/null | awk '{print $2}')

	# Track steps needed for fallback
        start_steps=0
        end_steps=0
	
	# If missing, try with multi-step context expansion
        if [[ -z "$ref_startpos" ]]; then
            result=$(find_ref_pos "$start_node" "$xg" "$ref_path_name")
            ref_startpos=$(echo "$result" | cut -d',' -f1)
            start_steps=$(echo "$result"  | cut -d',' -f2)
        fi

        if [[ -z "$ref_endpos" ]]; then
            result=$(find_ref_pos "$end_node" "$xg" "$ref_path_name")
            ref_endpos=$(echo "$result" | cut -d',' -f1)
            end_steps=$(echo "$result"  | cut -d',' -f2)
        fi

        echo -e "${scaffold}\t${id}\t${start_node}\t${ref_startpos}\t${start_steps}\t${end_node}\t${ref_endpos}\t${end_steps}" \
            >> ref_inversion_coordinates.tsv

    done < <(zcat "$vcf" | grep -v "^#")

done

echo "Done. Output in ref_inversion_coordinates.tsv"
