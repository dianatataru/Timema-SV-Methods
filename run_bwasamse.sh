#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/bwamem_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/bwamem_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=bwamem
#SBATCH --qos gompert-grn
#SBATCH --array=1-238   # Job array when n is number of unique samples

### LOAD MODULES ###
#For this step, bwa and needed
module load bwa
module load samtools

echo "Start Job"
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

### ASSIGN VARIABLES  ###
P=$(find /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/REF_data -type l | sort | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
SAMPLE=$(basename $P | cut -d "." -f 1)
pangenome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.sv.gfa.fa.gz"
genome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked"
echo "P=$P"
echo "SAMPLE=$SAMPLE"
echo "pangenome=$pangenome"
echo "genome=$genome"

### SET TMPDIR ###
WORKDIR="/scratch/general/nfs1/u6071015/timemaGBS/"
cd "$WORKDIR"

### MAPPING ###
echo "Mapping ${SAMPLE}"

#replicated Science paper
bwa aln -n 4 -k 2 -l 20 -q 10 "$genome" "$P" > "${SAMPLE}_aligned.sai"
bwa samse "$genome" "${SAMPLE}_aligned.sai" "$P" | \
    samtools view -bS -q 1 - | \
    samtools sort - > "${SAMPLE}.sorted.unique.bam"
samtools index "${SAMPLE}.sorted.unique.bam"

echo "Mapping complete for ${SAMPLE}"

