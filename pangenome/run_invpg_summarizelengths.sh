#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 6
#SBATCH --account=gompert
#SBATCH --partition=notchpeak
#SBATCH --job-name=summarizelengths
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/INVPG-annot/INVPG_annot/HWY154_REF_4119Hap2/summarizelengthsgenos%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/INVPG-annot/INVPG_annot/HWY154_REF_4119Hap2/summarizelengthsgenos-%j.out

module load bcftools

vcf="invpg_HWY154_REF_4119Hap2.vcf"

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/INVANNOT\t%INFO/AC\t%INFO/NS\t%INFO/AN\n' \
	"$vcf" > invpg_HWY154_REF_4119Hap2_all_v2.tsv

awk 'BEGIN{OFS="\t";
    # Print header with REF/ALT replaced by bp count columns
    print "CHROM","POS","ID","REF_bp","ALT_bp","INVANNOT","AC","NS","AN"
}
{
    # REF_bp: length of REF allele (col 4)
    ref_bp = length($4)

    # ALT_bp: comma-separated lengths for each ALT allele (col 5)
    n = split($5, alts, ",")
    alt_bp = ""
    for (i=1; i<=n; i++) {
        alt_bp = alt_bp (i>1 ? "," : "") length(alts[i])
    }

    # Print all fields, replacing REF and ALT with their bp counts
    print $1, $2, $3, ref_bp, alt_bp, $6, $7, $8, $9
}
' invpg_HWY154_REF_4119Hap2_all_v2.tsv > invpg_HWY154_REF_4119Hap2_sumbp_v2.tsv
