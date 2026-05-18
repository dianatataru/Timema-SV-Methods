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
echo -e "scaffold\tvariant_id\tgenome\tstart_node\tstart_pos\tend_node\tend_pos\tgraph_seq_len" > inversion_coordinates.tsv

# Loop over each scaffold's inversion VCF
for vcf in *_pantree_inversions_only.vcf.gz; do

    # Extract scaffold name/number from filename
    scaffold=$(basename "$vcf" .vcf.gz | sed 's/_pantree_inversions_only//')
    
    # Path to corresponding xg index, vg file, and display paths file
    xg="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/${scaffold}.xg"
    vg_file="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/${scaffold}.vg"
    paths_file="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/display_paths_${scaffold}.txt"
    
    # Check vg file exists
    if [[ ! -f "$vg_file" ]]; then
        echo "WARNING: $vg_file not found, skipping $scaffold" >&2
        continue
    fi

    # Build xg index if it doesn't exist
    if [[ ! -f "$xg" ]]; then
        echo "INFO: $xg not found, building from $vg_file..." >&2
        vg index -x "$xg" "$vg_file"
        if [[ $? -ne 0 ]]; then
            echo "ERROR: failed to build $xg, skipping $scaffold" >&2
            continue
        fi
        echo "INFO: $xg built successfully" >&2
    else
        echo "INFO: $xg already exists, skipping index build" >&2
    fi

    # Check display paths file exists
    if [[ ! -f "$paths_file" ]]; then
        echo "WARNING: $paths_file not found, skipping $scaffold" >&2
        continue
    fi

    # Loop over each inversion in the VCF (skip header lines, decompress on the fly)
    while IFS=$'\t' read -r chrom pos id ref alt qual filter info format samples; do
        
        # Extract start and end nodes from the ID field e.g. >2993340<2993702
        start_node=$(echo "$id" | grep -oP '(?<=[><])\d+' | head -1)
        end_node=$(echo "$id"   | grep -oP '(?<=[><])\d+' | tail -1)
        
	#extract length of the inversion
	graph_seq_len=$(vg find -x "$xg" -n "$start_node" -n "$end_node" -c 0 2>/dev/null | vg view -j - 2>/dev/null \
    		| python3 -c "import json,sys; g=json.load(sys.stdin); print(sum(len(n['sequence']) for n in g.get('node',[])))" 2>/dev/null)
	[[ -z "$graph_seq_len" ]] && graph_seq_len="."

        # Loop over each genome path
        while IFS= read -r genome_path; do
            
            # Skip empty lines
            [[ -z "$genome_path" ]] && continue
            
            # Query start node
            start_result=$(vg find -x "$xg" -n "$start_node" -P "$genome_path" 2>/dev/null)
            start_pos=$(echo "$start_result" | awk '{print $2}')
            
            # Query end node
            end_result=$(vg find -x "$xg" -n "$end_node" -P "$genome_path" 2>/dev/null)
            end_pos=$(echo "$end_result" | awk '{print $2}')
            
            # If either position is empty, mark as missing
            [[ -z "$start_pos" ]] && start_pos="."
            [[ -z "$end_pos" ]]   && end_pos="."
            
            # Write to output table
            echo -e "${scaffold}\t${id}\t${genome_path}\t${start_node}\t${start_pos}\t${end_node}\t${end_pos}\t${graph_seq_len}" \
                >> inversion_coordinates.tsv

        done < "$paths_file"

    done < <(zcat "$vcf" | grep -v "^#")

done

echo "Done. Output in inversion_coordinates.tsv"
