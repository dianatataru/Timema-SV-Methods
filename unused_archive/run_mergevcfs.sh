#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=mergevcf
#SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/merge-%A_%a.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/merge-%A_%a.out
#SBATCH --mem=200G

### Load Modules ###
module load bcftools
module load plink/2.0

#cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/vg_vcf

### merge vcfs ###
#bcftools merge *.vcf.gz -O z -o 2013FHA_merged.vcf.gz
#bcftools index 2013FHA_merged.vcf.gz


### check missingness ###
#plink2 --vcf 2013FHA_merged.vcf.gz --allow-extra-chr --missing --out 2013FHA_merged_missingness_report
#plink2 --vcf 2013FHA_merged.vcf.gz --allow-extra-chr --freq --out 2013FHA_merged_freq


# merge pantree output
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree

bcftools merge *_inversions_only.vcf.gz -O z -o pantree_inversions_only_merged.vcf.gz 
bcftools index pantree_inversions_only_merged.vcf.gz
