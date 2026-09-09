#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=summarizepantree
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizepantree-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizepantree-%j.out

module load miniforge3
cd /uufs/chpc.utah.edu/common/home/u6071015/software/pantree
source .venv/bin/activate
export PYTHONPATH="/uufs/chpc.utah.edu/common/home/u6071015/software/pantree:${PYTHONPATH}"

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree
SCAFF="Scaffold_13__3_contigs__length_82050896"

#mkdir summary

#summarize SVs (edited from pantree manuscript, puts output in /summary subdir of working directory)
python pantree_summary.py --vcf /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/${SCAFF}_pantree.vcf.gz --chrom ${SCAFF}


#subset vcf to only inversions
zcat ${SCAFF}_pantree.vcf.gz \
| awk '
  /^#/ { print; next }
  $8 ~ /(^|;)VT=INV(;|$)/
' \
| gzip > ${SCAFF}_pantree_inversions_only.vcf.gz
