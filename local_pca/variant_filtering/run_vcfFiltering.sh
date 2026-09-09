#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/filter_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/filter_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=filter
#SBATCH --qos gompert-grn

### LOAD MODULES ###
module load R
module load bcftools

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf/

#plot stats on raw vcf
bcftools query \
       -f '%CHROM\t%POS\t%QUAL\t%INFO/MQBZ\t%INFO/SCBZ\t%INFO/RPBZ\t%INFO/DP\n' \
      FHA_all_oneref.vcf > FHA_all_oneref_site_metrics.txt
Rscript plot_histograms.R

#filtered the SNP set for coverage, missing data, and various tests of bias 
perl vcfFilter.pl REF_all_oneref.vcf

#extracted the read depth per SNP and individual from the filtered vcf files
bcftools query -f '[%DP\t]\n' filtered2x_REF_all_oneref.vcf | sed 's/\t$//' > depth_oneref.txt

#compute depth per individual and SNP to identify SNPs and individuals to drop
Rscript CovFilt.R

#Drop those individuals
perl filterSomeMore.pl --keep KeepSNPs.txt --keepinds KeepInds.txt --out morefilter_2x_REF_all_oneref.vcf filtered2x_REF_all_oneref.vcf

#convert to genotype likelhood file
perl vcf2gl.pl 0.0 morefilter_2x_REF_all_oneref.vcf

#obtain maximum likelihood of AF using expectation-maximization algorithm written by Zach, estpEM
#before running, change first line to have #ind #loci instead of 0 0
./estpEM -i REF_all.gl -o REF_all_estpEM.txt -e 0.001 -m 50 -h 1
