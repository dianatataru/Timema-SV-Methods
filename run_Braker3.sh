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
busco -i /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked -m geno -o busco_genome_out -l insecta_odb10

## amino acids
busco -i braker.aa -m prot -o busco_aa_out -l insecta_odb10

## Augustus amino acids
cd Augustus #had to add this for it to find the input files
busco -i augustus.hints.aa -m prot -o busco_augustus_aa_out -l insecta_odb10
