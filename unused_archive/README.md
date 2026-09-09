# Unused Archive

These are scripts that were not used in the submitted manuscript, mostly having to do with alignment of genotyping-by-sequencing data to the pangenome. 
Additionally, some other analyses that were not included:

## GBS Data Alignment and Variant Calling from Pangenome with VG

*This is not in the paper*
We use the vg-giraffe-pack-call workflow to align and call SVs from the GSH2-8haplotype pangenome(https://link.springer.com/article/10.1186/s13059-020-1941-7#Sec12).

Input graph is the filtered pangenome output form MinigraphCactus that has been vg autoindexed to create the .dist, .shortread.zipcodes, and .shortread.withzip.min indexes. Snarls have been identified using vg snarls. This filter graph is made by removing nodes covered by fewer than 2 haplotypes from the clip graph. The distance index (.dist) is a memory-mapped file. As of vg version 1.48.0, the file will be opened in read+write mode by default. This can cause issues in HPC clusters and other distributed environments, where multiple computers try to access the same distance index file. To avoid this, make the file read-only or use a local copy of the file (chmod 444 HWY154_REF_4119Hap2.d2.dist). Joint variant calling is not possible with vg call. Here, I make a separate vcf for each sample and then will join them all with bcftools merge. I think it is possible to take the output .gam from giraffe and surject it to make bams, and then use GATK gvcf, although it is possible this would return nonsense (https://github.com/vgteam/vg/issues/1416). 

Output from cactus minigraph- generally the clip graph for everything except --giraffe which defaults to the filter graph, and anything odgi-related which defaults to full. 

Make soft linked GBS files:
```
ln -s /uufs/chpc.utah.edu/common/home/gompert-group3/data/sheffield/timema/2013fha_gwas/02_ids_reads/cristinae/2013*bz2 data/
```
now run alignment and variant calling:

```
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=mapGBS
#SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/mapGBS-%A_%a.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/mapGBS-%A_%a.out
#SBATCH --mem=200G
#SBATCH --array=0-1

module load cactus/3.0.1
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS

FILES=(data/*.fq.bz2)
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
id=$(basename "$FILE" .fq.bz2)
echo ID=${id}

PANGENOME_PATH="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2"
PANGENOME="HWY154_REF_4119Hap2.d2"

# Graph alignment (do I need -x graph.xg?
vg giraffe \
  -Z ${PANGENOME_PATH}/${PANGENOME}.gbz \
  -d ${PANGENOME_PATH}/${PANGENOME}.dist \
  -z ${PANGENOME_PATH}/${PANGENOME}.shortread.zipcodes \
  -m ${PANGENOME_PATH}/${PANGENOME}.shortread.withzip.min \
  -f <(bzcat data/${id}.fq.bz2) \
  -t 24 \
  > vg_intermediate/${id}.gam

vg stats -a vg_intermediate/${id}.gam > vg_stats/${id}.stats.txt

# Coverage packing
vg pack \
  -x ${PANGENOME_PATH}/${PANGENOME}.gbz \
  -g vg_intermediate/${id}.gam \
  -Q 5 \
  -t 20 \
  -o vg_intermediate/${id}.pack

# Variant calling
vg call \
  ${PANGENOME_PATH}/${PANGENOME}.gbz \
  -r ${PANGENOME_PATH}/HWY154_REF_4119Hap2.snarls \
  -k vg_intermediate/${id}.pack \
  -a -A --progress\
  -t 20 -z -c 50 -C 10000000\
  -s ${id} \
  > vg_vcf/${id}.vcf

```

output vcf:
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_4__1_contigs__length_97222829,length=97222830>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_9__2_contigs__length_79556474,length=79556476>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_1__1_contigs__length_160647932,length=160647933>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_6__1_contigs__length_78844258,length=78844259>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_13__3_contigs__length_82050896,length=82050899>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_3__2_contigs__length_137956696,length=137956698>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_2__1_contigs__length_157594471,length=157594472>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_12__1_contigs__length_47609450,length=47609451>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_11__2_contigs__length_80009992,length=80009994>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_8__1_contigs__length_71271319,length=71271320>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_7__1_contigs__length_75018798,length=75018799>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_5__1_contigs__length_83128659,length=83128660>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_10__2_contigs__length_75648701,length=75648703>
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=MAD,Number=1,Type=Integer,Description="Minimum site allele depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled probability estimate of the called genotype">
##FORMAT=<ID=GP,Number=1,Type=Float,Description="Genotype Probability, the log-scaled posterior probability of the called genotype">
##FORMAT=<ID=XD,Number=1,Type=Float,Description="eXpected Depth, background coverage as used for the Poisson model">
##FILTER=<ID=lowad,Description="Variant does not meet minimum allele read support threshold of 1">
##FILTER=<ID=lowdepth,Description="Variant has read depth less than 4">
##SAMPLE=<ID=2013FHA001>
##contig=<ID=Scaffold_10__2_contigs__length_75648701>
##contig=<ID=Scaffold_11__2_contigs__length_80009992>
##contig=<ID=Scaffold_12__1_contigs__length_47609450>
##contig=<ID=Scaffold_13__3_contigs__length_82050896>
##contig=<ID=Scaffold_1__1_contigs__length_160647932>
##contig=<ID=Scaffold_2__1_contigs__length_157594471>
##contig=<ID=Scaffold_3__2_contigs__length_137956696>
##contig=<ID=Scaffold_4__1_contigs__length_97222829>
##contig=<ID=Scaffold_5__1_contigs__length_83128659>
##contig=<ID=Scaffold_6__1_contigs__length_78844258>
##contig=<ID=Scaffold_7__1_contigs__length_75018798>
##contig=<ID=Scaffold_8__1_contigs__length_71271319>
##contig=<ID=Scaffold_9__2_contigs__length_79556474>
##bcftools_viewVersion=1.23-3-g34a49760+htslib-1.23-9-gacc28ac1
##bcftools_viewCommand=view -h 2013FHA001.vcf.gz; Date=Wed Mar  4 19:32:45 2026
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2013FHA001
```
bcftools +counts for same vcf:
```
bcftools +counts 2013FHA001.vcf.gz
[W::bcf_hdr_check_sanity] GP should be declared as Number=G
Number of samples: 1
Number of SNPs:    29859
Number of INDELs:  239728
Number of MNPs:    36580
Number of others:  78012
Number of sites:   413902
```

then can use bcftools merge to merge all vcfs? Do I need to rpoject this back onto snarls somehow to identify inversions? when I look at the snarls.json, there are 45265393 lines. Of those, there are 11779457 instances of the word "backward" which I think is related to the inversions because when I search one of the nodes of the inversion (9179874) it pops up in a set of parenthesis including backwards=true.
to merge vcfs:
```
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=mergevcf
#SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/merge-%A_%a.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/merge-%A_%a.out
#SBATCH --mem=200G

### load modules ###
module load bcftools
module load plink/2.0 

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/vg_vcf

### merge and index vcfs ###
bcftools merge *.vcf.gz -O z -threads 12 -o 2013FHA_merged.vcf.gz
bcftools index 2013FHA_merged.vcf.gz

### check missingness ###
plink --vcf 2013FHA_merged.vcf.gz --missing --out 2013FHA_merged_missingness_report
```
## Genome Annotation and GENESPACE visualization

*This was not actually used in the paper but I am keeping the section for later reference*
  
We can use the genespace visualization to validate the inversions and translocations found.
Copying over the braker3 annotations from the Science Paper:

```
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation

cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h1/braker/braker.aa t_crist_refug_green_h1.fa #REDOWNLOAD THIS
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_stripe_h1/brakerV1/braker.aa t_crist_refug_stripe_h1.fa #not found, changed dir name
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h2/braker/braker.aa t_crist_h154_green_h2.fa #not found, reran braker but OLD
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h1/braker/braker.aa t_crist_h154_stripe_h1.fa
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h1/brakerV1/braker.aa t_crist_h154_green_h1.fa #not found, changed dir name
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h2/braker/braker.aa t_crist_refug_green_h2.fa
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_stripe_h2/braker/braker.aa t_crist_refug_stripe_h2.fa
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h2/brakerV1/braker.aa t_crist_h154_stripe_h2.fa #not found, changed dir name


## fix format
perl -p -i -e 's/\.t[0-9]//' *fa

grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h1/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_refug_green_h1.bed
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h2/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_refug_green_h2.bed
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_stripe_h1/brakerV1/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_refug_stripe_h1.bed 
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h2/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_h154_green_h2.bed #not made
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h1/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_h154_stripe_h1.bed 
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h1/brakerV1/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_h154_green_h1.bed 
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_stripe_h2/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_refug_stripe_h2.bed 
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h2/brakerV1/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_h154_stripe_h2.bed 
```

One of the annotations is missing so need to run braker on it

```
#!/bin/bash 
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=braker
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/braker-%j.err
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/braker-%j.out

source ~/.bashrc

ml braker/3.0.8
ml busco

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/t_crist_hyw154_green_h2

## run braker

braker.pl --genome=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked \
	--prot_seq=/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/proteins.fasta \
	--rnaseq_sets_ids=clean_tcr135.17_0003_R,clean_tcr137.17_0006_R,clean_tcr139.17_0012_R,clean_tcr140.17_0015_R,clean_tcr141.17_0019_R,clean_tcr142.17_0043_R,clean_tcr143.17_0045_R,clean_tcr144.17_0049_R,clean_tcr145.17_0051_R,clean_tcr146.17_0057_R,clean_tcr148.17_0062_R,clean_tcr149.17_0065_R,clean_tcr150.17_0067_R,clean_tcr151.17_0070_R,clean_tcr152.17_0074_R,clean_tcr173.17_0075_R,clean_tcr174.17_0081_R,clean_tcr175.17_0082_R \
	--rnaseq_sets_dirs=/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/rna_seq_for_annotations \
	--AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts \
	--AUGUSTUS_CONFIG_PATH=/uufs/chpc.utah.edu/common/home/u6071015/augustus/config \
	--threads=48 --gff3

## run busco, genome and aa
cd braker
## genome
#busco -i /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked -m geno -o busco_genome_out -l insecta_odb10

## amino acids
busco -i braker.aa -m prot -o busco_aa_out -l insecta_odb10

## Augustus amino acids
cd Augustus #had to add this for it to find the input files
busco -i augustus.hints.aa -m prot -o busco_augustus_aa_out -l insecta_odb10
```

then prepare it like the rest of the files:
```
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation
cp /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/t_crist_hyw154_green_h2/braker/braker.aa t_crist_h154_green_h2.fa
perl -p -i -e 's/\.t[0-9]//' t_crist_h154_green_h2.fa
grep "gene" /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/t_crist_h154_green_h2/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_hyw154_green_h2.bed

#there are peptides missing from the bed genes, which means I need to filter the .fa files to match
module load seqkit/2.8.2
SAMPLE="t_crist_h154_green_h2"
cut -f4 ${SAMPLE}.bed | sort -u > bed.ids
seqkit grep -f bed.ids ${SAMPLE}.fa -o ${SAMPLE}_filtered.fa
cp ${SAMPLE}_filtered.fa /scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/peptide/${SAMPLE}.fa
#for samples where the filtered matched unfiltered, I deleted the filtered file.
```

then we can run GENESPACE:

first make folders and move files to those folders:

```
cd /scratch/general/nfs1/u6071015
mkdir GENESPACE_TIMEMA
cd GENESPACE_TIMEMA/
mkdir peptide
mkdir bed
cp /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/*filtered.fa peptide/
cp /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/*bed bed/

```
then run the program
```
#!/bin/bash 
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos=gompert-grn
#SBATCH --job-name=GENESPACE
#SBATCH --error=/scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/GENESPACE-%j.err
#SBATCH --output=/scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/GENESPACE-%j.out

#load modules
module load orthofinder
module load R

cd /scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/

echo "start GENESPACE"

Rscript genespace_TIMEMA.R

echo "GENESPACE done"
```

and this is genespace_TIMEMA.R:

```
#create personal library to write packages to
#dir.create("~/R/x86_64-pc-linux-gnu-library/4.4", recursive = TRUE, showWarnings = FALSE)
#set library paths
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")

#only run devtools download for the first time
#if (!requireNamespace("devtools", quietly = TRUE))
#    install.packages("devtools")
#devtools::install_github("jtlovell/GENESPACE")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("Biostrings", "rtracklayer"))

library(GENESPACE)

#set working directory
wd<-"/scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/"
path2mcscanx<-"~/bin/MCScanX/"

# initalize the run and QC the inputs
gpar<-init_genespace(wd=wd,path2mcscanx=path2mcscanx)

# need to set this
gpar$shellCalls$orthofinder<-"orthofinder"

# accomplish the run
out <- run_genespace(gpar, overwrite = T)

# plot
roi<-data.frame(
		genome=c("t_crist_h154_green_h1","t_crist_h154_green_h2",
			 "t_crist_h154_stripe_h1","t_crist_h154_stripe_h2",
		"t_crist_refug_green_h1","t_crist_refug_green_h2",
    "t_crist_refug_stripe_h1","t_crist_refug_stripe_h2"),
		start=c(0,0,0,0,0,0,0,0,0,0),end=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf))

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))

pdf("syn8way.pdf",width=9,height=5.6)
ripd <- plot_riparian(
	gsParam = out,
	palette = customPal,
        #highlightBed = roi,
	braidAlpha = .3,
	useOrder=TRUE,
	chrFill = "lightgrey",
    addThemes = ggthemes,
	useRegions = FALSE,
	    #invertTheseChrs = invchr,
  	refGenome = "t_crist_h154_stripe_h2",
	genomeIDs = c("t_crist_h154_green_h1","t_crist_h154_green_h2",
			 "t_crist_h154_stripe_h1","t_crist_h154_stripe_h2",
		"t_crist_refug_green_h1","t_crist_refug_green_h2",
    "t_crist_refug_stripe_h1","t_crist_refug_stripe_h2"),
	backgroundColor = NULL)
dev.off()
```

## SV calling from comparative alignment with vg
Ultimately this did not really work but I could see this being useful down the line.
halSummarizeMutations also doesn't come up with identifiers for mutations, so what I really need is a SV caller that develops position/reference-specific (?) identifiers so that I can then tell how many are unique across the pairwise comparisons. The trick is to do this without accidentally just creating another pangenome... it seems like vg might have a way of calling SVs which is rigorous compared to other SV callers (Hickey et al. Genome Biology (2020) https://doi.org/10.1186/s13059-020-1941-7).

New Jay paper does the following with vg deconstruct vcf output:
- ran vcfbub to keep only top-level variant sites (snarls) less than 100 kb in size
- used vcfwave to realign REF and ALT alleles to split nested alleles to separate entries and identify inversions >1kb
- combined vcf files with bcftools concat, added in missing sample coolumns with bcftools query, and used bcftools fixploidy to set allele number for every site and bcftools fill tags to add AF and AC for each each site. Also used bcftools norm to split multiallelic to biallelic
  
For the following SV calling, GSH2 is the REF for all. 
```
salloc --time=06:00:00 --ntasks 1 --nodes=1 --account=gompert --partition=gompert-grn --qos gompert-grn --mem=300G
module load cactus/3.0.1
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus

hal2vg cactusStripe_TcrGSH2_TcrGUSH2_DTv2.hal --hdf5InMemory --chop 32 --progress > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg
vg index cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg -L
vg snarls cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls
vg view -j -R cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls.json
wc -l cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls.json
#14234556 cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls.json

vg deconstruct cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg -P TcrGSH2 -e -a > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vcf
bgzip cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vcf
tabix cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vcf.gz

module load bcftools
# subset to SVs greater than 0 bp, throws up error about GT but just rem0ves it
bcftools view -i 'strlen(REF)>50 || strlen(ALT)>50' \
    cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vcf.gz -G -Oz -o cactusStripe_TcrGSH2_TcrGUSH2_min50bp.vcf.gz
tabix -p vcf cactusStripe_TcrGSH2_TcrGUSH2_min50bp.vcf.gz
bcftool stats cactusStripe_TcrGSH2_TcrGUSH2_min50bp.vcf.gz # number of records is 819, I think this is missing a lot of stuff

#truncated due to this VCF parse error:
#Couldn't read GT data: value not a number or '.' at TcrGSH2#0#Scaffold_10__2_contigs__length_75648701:281081 and #TcrGSH2#0#Scaffold_10__2_contigs__length_75648701:2810840
#investigate just that position to see what is wrong
vg chunk \
  -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg \
  -p "TcrGSH2#0#Scaffold_10__2_contigs__length_75648701:2810800-2810900" \
  --snarls cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls \
  -g > 2810800-2810900.vg

vg view -d 2810800-2810900.vg > 2810800-2810900.dot
dot -Tpdf 2810800-2810900.dot > 2810800-2810900.pdf

vg paths -v 2810800-2810900.vg -Q TcrGSH2 -L
TcrGSH2#0#Scaffold_10__2_contigs__length_75648701[2810781]
TcrGSH2#0#Scaffold_13__3_contigs__length_82050896[3492648]
TcrGSH2#0#Scaffold_10__2_contigs__length_75648701[6014359]

# this is a translocation, which deconstruct can't handle. I actually want to be using vg call.
# First I need to make associated indexes and a gam for each fasta

vg index -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg \
         -g cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gcsa \
         -L -j cactusStripe_TcrGSH2_TcrGUSH2_DTv2.dist \
         cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg

vg minimizer cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg \
  -d cactusStripe_TcrGSH2_TcrGUSH2_DTv2.dist \
  -o cactusStripe_TcrGSH2_TcrGUSH2_DTv2.min

vg giraffe \
  -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg \
  -g cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gcsa \
  -m cactusStripe_TcrGSH2_TcrGUSH2_DTv2.min \
  -d cactusStripe_TcrGSH2_TcrGUSH2_DTv2.dist \
  -f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked \
  > TcrGSH2.gaf

vg giraffe \
  -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg \
  -g cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gcsa \
  -m cactusStripe_TcrGSH2_TcrGUSH2_DTv2.min \
  -d cactusStripe_TcrGSH2_TcrGUSH2_DTv2.dist \
  -f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked \
  > TcrGUSH2.gaf

vg prune -r -p -t 2 cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.pruned.vg
vg index cactusStripe_TcrGSH2_TcrGUSH2_DTv2.pruned.vg -L -j cactusStripe_TcrGSH2_TcrGUSH2_DTv2.pruned.dist 

vg autoindex --workflow giraffe -g cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gfa \
	-p cactusStripe_TcrGSH2_TcrGUSH2_DTv2 \
	-G cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gbz \
	--threads 2 --target-mem 5G --verbosity 2 -T temp/

# Constructing distance index for Giraffe.
#Killed
```
or as sbatch script 
```
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus
#SBATCH --qos gompert-grn
#SBATCH -e cactus-%j.err
#SBATCH -o cactus-%j.out
#SBATCH --mem=300G

module load cactus/3.0.1
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus

export TMPDIR="/scratch/general/nfs1/u6071015/vg_tmp"

PAIR="cactusStripe_TcrGSH2_TcrGUSH2_DTv2"
GENOME1="t_crist_hwy154_cen4119_hap2.fasta.masked"
GENOME2="t_crist_hwy154_cen4280_hap2.fasta.masked"
SAMP1="TcrGSH2"
SAMP2="TcrGUSH2"

#hal2vg ${PAIR}.hal --hdf5InMemory --chop 32 --progress > ${PAIR}.vg

#vg index ${PAIR}.vg -x ${PAIR}.xg -L

#vg snarls ${PAIR}.xg > ${PAIR}.snarls

#vg convert -f ${PAIR}.vg > ${PAIR}.gfa

#vg gbwt --num-jobs 16 --gbz-format -g ${PAIR}.gbz -G ${PAIR}.gfa

vg prune -r -p -t 16 ${PAIR}.vg > ${PAIR}.pruned.vg

vg index ${PAIR}.vg \
         -L -j ${PAIR}.dist 

vg giraffe -b hifi -Z ${PAIR}.gbz \
	-f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/${GENOME1} \
	-p > ${SAMP1}.gam

vg giraffe -b hifi -Z ${PAIR}.gbz \
	-f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/${GENOME2} \
	-p > ${SAMP2}.gam

vg pack ${PAIR}.vg \
        -g ${SAMP1}.gam \
        -g ${SAMP2}.gam \
        -o ${PAIR}.pack

vg call -A -c 50 -r ${PAIR}.snarls \
	--threads 6 -S ${SAMP1} \
	-k ${PAIR}.pack \
	${PAIR}.vg > ${PAIR}.vcf.gz

```
vg prune worked with a temp dir added, but vg index still oom killing after:

```
[vg prune] Original graph cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg: 75689771 nodes, 90722944 edges
[vg prune] Built a temporary XG index
[vg prune] Removed all paths
[vg prune] Pruned complex regions: 75689771 nodes, 81196671 edges
[vg prune] Removed small subgraphs: 66856969 nodes, 77460277 edges
Restored graph: 75689771 nodes
[vg prune] Serialized the graph: 75689771 nodes, 90722944 edges
INFO:    gocryptfs not found, will not be able to use gocryptfs
/uufs/chpc.utah.edu/sys/installdir/lmod/8.6-r8/init/bash: line 82: 1265292 Killed                  apptainer exec --nv /uufs/chpc.utah.edu/sys/installdir/cactus/3.1.4/cactus-3.1.4.sif vg $@
slurmstepd: error: Detected 1 oom_kill event in StepId=814509.batch. Some of the step tasks have been OOM Killed.
```

## Trying a different way of calling SV from pangenome 

Using the program INVPG_annot (https://github.com/SandraLouise/INVPG_annot)

```
#install in /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/INVPG-annot
git clone https://github.com/SandraLouise/INVPG_annot.git
cd INVPG_annot
pip install -r requirements.txt --upgrade
python -m pip install . --quiet
```

```
salloc --time=10:00:00 --ntasks 12 --nodes=1 --account=gompert-np --partition=gompert-np
salloc --time=10:00:00 --ntasks 1 --nodes=1 --account=gompert --partition=gompert-grn --qos=gompert-grn

#usage: invpg [-h] [-v INPUT_VCF_FILE] [-g INPUT_GFA_FILE] [-o OUTPUT_PREFIX] [-d DIV_PERCENTAGE] [-m MINCOV] [-k] [-t THREADS] [-
#test
cd test-dir
invpg -v test_bubbles.vcf -g test_graph.gfa -o test_annotation.vcf -m 0.5 -d 10
diff expected_annotation.vcf test_annotation.vcf

#softlink into HWY154_REF_4119Hap2/ folder
ln -s /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.vcf
ln -s /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.gfa

#run
invpg  -v HWY154_REF_4119Hap2.vcf -g HWY154_REF_4119Hap2.gfa -o invpg_HWY154_REF_4119Hap2.vcf -m 0.5 -d 10 -t 12
#Bubbles after filtering: 209656
#Inversion annotated bubbles: 403
#Results output in files invpg_HWY154_REF_4119Hap2.vcf.vcf and invpg_HWY154_REF_4119Hap2.vcf.stats

Total_bubbles   33224918
Large_bubbles   209656
Inversion_bubbles       403
Path-explicit   45
Alignment-rescued       519
```
then to extract inversion information from the vcf:

```
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=summarizelengths
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/INVPG-annot/INVPG_annot/HWY154_REF_4119Hap2/summarizelengthsgenos%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/INVPG-annot/INVPG_annot/HWY154_REF_4119Hap2/summarizelengthsgenos-%j.out

module load bcftools

vcf="invpg_HWY154_REF_4119Hap2.vcf"

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/INVANNOT\t%INFO/AC\t%INFO/NS\t%INFO/AN\n' invpg_HWY154_REF_4119Hap2.vcf > invpg_HWY154_REF_4119Hap2_all.tsv

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
' invpg_HWY154_REF_4119Hap2_all.tsv > invpg_HWY154_REF_4119Hap2_sumbp.tsv
```



