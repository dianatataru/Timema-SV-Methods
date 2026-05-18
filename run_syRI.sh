#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos=gompert-grn
#SBATCH --job-name=syri
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/syri-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/syri-%j.out

module load miniforge3
module load minimap2/2.24
conda activate syRI

#set paths and variables
QRY="TcrRGUS1"
cwd="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/syri" 
REFGENOME="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_hwy154_cen4119_hap2.fasta.masked.fused_RGUS1"
QRYGENOME="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_refug_cen4120_hap1.fasta.masked.reoriented"
OUT="syri_TcrGSH2_${QRY}"

cd ${cwd}
mkdir syri_TcrGSH2_${QRY}
cd syri_TcrGSH2_${QRY}

echo -e "${REFGENOME}\tHGS2\n${QRYGENOME}\tRGUS1" > genomes.txt

#perform whole genome alignment
minimap2 -ax asm5 --eqx ${REFGENOME} ${QRYGENOME} > ${OUT}.sam

#runSyRI, -k means keep intermediate files, -F is input format, .sam
syri -c ${OUT}.sam -r ${REFGENOME} -q ${QRYGENOME} -k -F S --nosnp 

#Plotting genomic structures predicted by SyRI
plotsr --sr syri.out \
	--genomes genomes.txt \
	-H 8 -W 5 \
	-o ${OUT}_plot.pdf
