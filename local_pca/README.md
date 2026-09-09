# Local PCA

## Variant Calling
This folder contains the scripts that takes input of raw genotyping-by-sequencing data and creates output of variants for local PCA analysis \
- ```run_bwasamse.sh``` map samples to the reference genome with bwa aln/samse \
- ```run_varcallbcftools.sh``` joint variant calling with bcftools concensus caller \
- ```run_vcfFiltering.sh``` variant filtering \
  - ```plot_histograms.R ``` plot raw stats\
  - ```vcfFilter.pl``` filtered the SNP set for coverage, missing data, and various tests of bias  \
  - ```CovFilt.R```  compute depth per individual and SNP to identify SNPs and individuals to drop \
  - ```filterSomeMore.pl```  Drop those individuals \
  - ```vcf2gl.pl``` convert vcf to genotype likelihood file \
- ```run_gl2genest.sh``` Bayesian estimate of genotype likehood

## Local PCA analysis
- ```run_lostruct.sh``` this batch script creates a positions file from the vcf which is needed to run the localpca_manyMDSaxes_withfiltering.R script \
- ```localpca_manyMDSaxes_withfiltering.R``` this script takes the cpntest.txt and positions.txt file, and runs the local PCA test, flags outliers,
  conducts filtering based off of length and clustering, and plots significant outliers in both per-chromosome MDS plots and PCAs of the outlier region.
