#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH --mem=200G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=odgi-viz
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/logs/odgi-viz-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/logs/odgi-viz-%j.out

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate odgi

WORKDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments"
cd $WORKDIR

SCAFFOLDS=(
  Scaffold_1__1_contigs__length_160647932
  Scaffold_2__1_contigs__length_157594471
  Scaffold_3__2_contigs__length_137956696
  Scaffold_4__1_contigs__length_97222829
  Scaffold_5__1_contigs__length_83128659
  Scaffold_6__1_contigs__length_78844258
  Scaffold_7__1_contigs__length_75018798
  Scaffold_8__1_contigs__length_71271319
  Scaffold_9__2_contigs__length_79556474
  Scaffold_10__2_contigs__length_75648701
  Scaffold_11__2_contigs__length_80009992
  Scaffold_13__3_contigs__length_82050896
)

for SCAFFOLD in "${SCAFFOLDS[@]}"; do
    echo "Processing ${SCAFFOLD}"

vg convert -f ${SCAFFOLD}.vg -W > ${SCAFFOLD}.gfa

odgi build -g ${SCAFFOLD}.gfa  -o ${SCAFFOLD}.og -O

odgi sort -i ${SCAFFOLD}.og --threads 20 -P -Y -C /scratch/general/nfs1/u6071015/odgi -o ${SCAFFOLD}_sortedPGSGD.og

#make a file that lists just the genome paths. the graph will be too big with all MINIGRAPH paths.
mapfile -t TCRIST_PATHS < <(
  odgi paths -L -i ${SCAFFOLD}_sortedPGSGD.og | grep 't_crist'
)
PATHFILE="display_paths_${SCAFFOLD}.txt"
printf "%s\n" "${TCRIST_PATHS[@]}" > ${PATHFILE}

odgi viz \
    -i ${SCAFFOLD}_sortedPGSGD.og \
    -o ${SCAFFOLD}_sortedPGSGD.png \
    -p ${PATHFILE}

done
