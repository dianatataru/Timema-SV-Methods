#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=summarizelengths
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizelengthsgenos%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizelengthsgenos-%j.out

module load bcftools

FIRST=1

for vcf in *_pantree_inversions_only.vcf.gz; do
    scaffold=$(basename "$vcf" .vcf.gz | sed 's/_pantree_inversions_only//')

    SAMPLES=$(bcftools query -l "$vcf" \
      | grep 't_crist' \
      | grep -v 'MINIGRAPH' \
      | tr '\n' ',' \
      | sed 's/,$//')

    if [[ -z "$SAMPLES" ]]; then
        echo "WARNING: no t_crist samples found in $vcf, skipping" >&2
        continue
    fi

    if [[ $FIRST -eq 1 ]]; then
        # Print header on first VCF only
        bcftools query -s "$SAMPLES" -H \
          -f '%CHROM\t%ID\t%INFO/RC\t%INFO/AC\t%INFO/TP\t%INFO/NIA\t%INFO/AN\t%INFO/NR[\t%GT:%CR:%CA]\n' \
          "$vcf"
        FIRST=0
    else
        # No header for subsequent VCFs
        bcftools query -s "$SAMPLES" \
          -f '%CHROM\t%ID\t%INFO/RC\t%INFO/AC\t%INFO/TP\t%INFO/NIA\t%INFO/AN\t%INFO/NR[\t%GT:%CR:%CA]\n' \
          "$vcf"
    fi
done \
| awk 'BEGIN{OFS="\t"}
    /^#/ {
        sub(/NR/, "NR\tNR_count")
        print; next
    }
    {
        nr_count = ($8 != ".") ? 1 : 0
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", \
            $1, $2, $3, $4, $5, $6, $7, $8, nr_count
        for (i=9; i<=NF; i++) printf "\t%s", $i
        printf "\n"
    }' > all_scaffolds_inversions_tcrist_genotypes_allNR.tsv

awk 'BEGIN{OFS="\t"}
  NR==1 {
    # Fix header: remove NR col (8), rename NR_count col (9)
    for (i=1; i<=NF; i++) {
      if (i==8) continue
      if (i==9) printf "NR_bp"
      else printf "%s", $i
      if (i!=NF) printf "\t"
    }
    printf "\n"
    next
  }
  {
    nr_val = $8
    nr_bp = (nr_val == ".") ? 0 : length(nr_val)
    for (i=1; i<=NF; i++) {
      if (i==8) continue
      if (i==9) printf "%s", nr_bp
      else printf "%s", $i
      if (i!=NF) printf "\t"
    }
    printf "\n"
  }
' all_scaffolds_inversions_tcrist_genotypes_allNR.tsv > all_scaffolds_inversions_tcrist_genotypes_NRbp.tsv
