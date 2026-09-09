#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH --mem=200G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=odgi-draw
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/logs/odgi-draw-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/logs/odgi-draw-%j.out

module load cactus/3.0.1

cactus-graphmap-join /scratch/general/nfs1/u6071015/cactusNp/timema/timemadrawJSscaff9 \
   --vg /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/Scaffold_9__2_contigs__length_79556474.vg \
   --reference Hap2_t_crist_hwy154_cen4119 \
   --outDir /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/chroms \
   --outName HWY154_REF_4119Hap2_Scaff9 --draw 
