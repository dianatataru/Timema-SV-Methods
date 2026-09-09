#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/varcall_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/varcall_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=varcall
#SBATCH --qos gompert-grn
#SBATCH --mem=400G

### LOAD MODULES ###
module load samtools/1.16
module load bcftools/1.16
#ran previously with bcftools/1.23

### ASSIGN VARIABLES ###
BAMDIR="/scratch/general/nfs1/u6071015/timemaGBS/"
BAM_FILES=($(find "$BAMDIR" -type f -name "*.sorted.unique.bam" | sort ))
pangenome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.sv.gfa.fa.gz"
genome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked"
WORKDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf"
THREADS=12
MERGED="${WORKDIR}/REF_all_oneref.unique.bam"
SORTED="${WORKDIR}/REF_all_oneref.sorted.unique.bam"
OUTVCF="REF_all_oneref.vcf"

### MERGE ALL BAMS FOR VARIANT CALLING ###
echo "Merge BAM files"
cd "$WORKDIR"

samtools merge -f -r -c -p -@ ${THREADS} "$MERGED" "${BAM_FILES[@]}"
samtools sort -@ 12 -o "$SORTED" "$MERGED"
samtools index "$SORTED"
samtools flagstat "$SORTED"

echo "BAM files merged"

### VARIANT CALLING ###
echo "start variant calling"
cd "$WORKDIR"

#filter for one ref
bcftools mpileup -Ou -d 100000 -a DP,AD,ADF,ADR -q 20 -Q 30 -f "$genome" "$SORTED" | bcftools call -v -c -p 0.01 -Ov -o "$OUTVCF"

#no filter for pangenome
bcftools mpileup -Ou -d 100000 -a DP,AD,ADF,ADR -f "$genome" "$SORTED" | bcftools call -v -c -Ov -o "$OUTVCF"

#same commands as science paper
bcftools mpileup -Ou -d 500 -a DP,AD,ADF,ADR -Q 30 -q 20 --skip-indels -f "$genome" "${BAM_FILES[@]}" | \
    bcftools call -v -c -p 0.01 -P 0.001 -Ov -o "$OUTVCF"

echo "finished variant calling"
