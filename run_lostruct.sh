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
#bcftools query -f '%CHROM\t%POS\n' morefilter_2x_REF_all_oneref.vcf | sed 's/.*|s//' > positions.txt

#old script
#Rscript ${SCRIPTDIR}/lostruct_oneref.R

#new script
#Rscript ${SCRIPTDIR}/localpca.R cpntest_FHA_all.txt FHA_all

#updated for many MDS axes
#Usage: Rscript ${SCRIPTDIR}/localpca_manyMDSaxes_v2.R <input_file> <output_prefix> <window_size_snps> <n_axes>
Rscript ${SCRIPTDIR}/localpca_manyMDSaxes_v3.R cpntest_FHA_all.txt FHA_all 100 10
