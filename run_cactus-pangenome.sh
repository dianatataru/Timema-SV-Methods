#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus-pangenome
#SBATCH --qos=gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus-pangenome-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus-pangenome-%j.out

module load cactus/3.0.1

cd /scratch/general/nfs1/u6071015/cactusNp/timema/

#full wrapper
cactus-pangenome timema8hapsJS \
  /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/HWY154_REF.txt \
  --outDir /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome \
  --outName HWY154_REF_4119Hap2 \
  --reference Hap2_t_crist_hwy154_cen4119.2 \
  --maxCores 6 \
  --vcf --giraffe --gfa --gbz --viz

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome

#Summarize Mutations
halSummarizeMutations HWY154_4119Hap2.full.hal
#halSummarizeMutations HWY154.full.hal
