#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos=gompert-grn
#SBATCH --job-name=gretl
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/gretl-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/gretl-%j.out

module load miniforge3
conda activate gretl_env

INPUT_DIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments"
WORK_DIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/gretl"

SCAFFS=(
#  Scaffold_1__1_contigs__length_160647932
#  Scaffold_2__1_contigs__length_157594471
#  Scaffold_3__2_contigs__length_137956696
#  Scaffold_4__1_contigs__length_97222829
#  Scaffold_5__1_contigs__length_83128659
#  Scaffold_6__1_contigs__length_78844258
#  Scaffold_7__1_contigs__length_75018798
#  Scaffold_8__1_contigs__length_71271319
#  Scaffold_9__2_contigs__length_79556474
#  Scaffold_10__2_contigs__length_75648701
#  Scaffold_11__2_contigs__length_80009992
  Scaffold_12__1_contigs__length_47609450 
#  Scaffold_13__3_contigs__length_82050896
)

for SCAFF in "${SCAFFS[@]}"; do
    echo "Processing ${SCAFF}"

echo "Graph-based and hybrid stats"
gretl stats -g ${INPUT_DIR}/${SCAFF}.gfa -o ${WORK_DIR}/gretl_stats_${SCAFF}.txt

echo "Path-based statistics"
gretl stats -g ${INPUT_DIR}/${SCAFF}.gfa -o ${WORK_DIR}/gretl_pathstats_${SCAFF}.txt -p

done

echo "job done"
