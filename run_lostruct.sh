#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/lostruct_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/lostruct_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=lostruct
#SBATCH --qos gompert-grn

module load bcftools
module load R

SCRIPTDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf"
WORKDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf/FHA_alignedtocen4119hap2"

cd ${WORKDIR}

#make positions file
#bcftools query -f '%CHROM\t%POS\n' morefilter_2x_REF_all_oneref.vcf | sed 's/.*|s//' > positions_raw.txt
#awk -F'\t' 'BEGIN {
#    map[1]=4; map[2]=3; map[3]=13; map[4]=8; map[5]=6;
#    map[6]=2; map[7]=10; map[8]=7; map[9]=9; map[10]=11;
#    map[11]=12; map[12]=5; map[13]=1
#}
#{
#    match($1, /Scaffold_([0-9]+)__/, arr)
#    n = arr[1]+0
#    $1 = (n in map) ? "Chr" map[n] : "UNMAPPED_" $1
#    print $1 "\t" $2
#}' positions_raw.txt > positions.txt


#Usage: Rscript ${SCRIPTDIR}/localpca_manyMDSaxes_v2.R <input_file> <output_prefix> <window_size_snps> <n_axes>
Rscript ${SCRIPTDIR}/localpca_manyMDSaxes_withfiltering.R cpntest_FHA_all.txt FHA_all 100 40
#Rscript ${SCRIPTDIR}/localpca_manyMDSaxes_withfiltering.R cpntest_REF_all.txt REF_all 100 40
