#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=200G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos=gompert-grn
#SBATCH --job-name=pantree
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/pantree-%j.err
##SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/pantree-%j.out

SCAFF="Scaffold_13__3_contigs__length_82050896"

module load miniforge3

#convert vg
#conda activate odgi
#cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments
#vg convert -f  ${SCAFF}.vg -W > ${SCAFF}.gfa
#conda deactivate

#run pantree
cd /uufs/chpc.utah.edu/common/home/u6071015/software/pantree
source .venv/bin/activate
export PYTHONPATH="/uufs/chpc.utah.edu/common/home/u6071015/software/pantree:${PYTHONPATH}"

uv run pantree /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/${SCAFF}.gfa \
	/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/${SCAFF}_pantree.vcf.gz \
	--ref-name Hap2_t_crist_hwy154_cen4119.2 \
	--log-path /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/${SCAFF}.log \
	--chr-id ${SCAFF} \
	--priority-samples t_crist_hwy154_cen4119.1,t_crist_hwy154_cen4280.1,t_crist_hwy154_cen4280.2,t_crist_refug_cen4122.1,t_crist_refug_cen4122.2,t_crist_refug_cen4120.1,t_crist_refug_cen4120.2

