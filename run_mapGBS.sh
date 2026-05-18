#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=mapGBS
#SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/mapGBS-%A_%a.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/mapGBS-%A_%a.out
#SBATCH --array=5-602
#SBATCH --mem=100G

module load cactus/3.0.1
module load bcftools
module load htslib

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS

FILES=(data/*.fq.bz2)
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
id=$(basename "$FILE" .fq.bz2)
echo ID=${id}

PANGENOME_PATH="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2"
PANGENOME="HWY154_REF_4119Hap2.d2"

# Graph alignment
#vg giraffe \
#  -Z ${PANGENOME_PATH}/${PANGENOME}.gbz \
#  -d ${PANGENOME_PATH}/${PANGENOME}.dist \
#  -z ${PANGENOME_PATH}/${PANGENOME}.shortread.zipcodes \
#  -m ${PANGENOME_PATH}/${PANGENOME}.shortread.withzip.min \
#  -f <(bzcat data/${id}.fq.bz2) \
#  -t 24 \
#  > vg_intermediate/${id}.gam

#get stats
#vg stats -a vg_intermediate/${id}.gam > vg_stats/${id}.stats.txt

# Coverage packing
#vg pack \
#  -x ${PANGENOME_PATH}/${PANGENOME}.gbz \
#  -g vg_intermediate/${id}.gam \
#  -Q 5 \
#  -t 20 \
#  -o vg_intermediate/${id}.pack

# Variant calling
vg call \
  ${PANGENOME_PATH}/${PANGENOME}.gbz \
  -r ${PANGENOME_PATH}/HWY154_REF_4119Hap2.snarls \
  -k vg_intermediate/${id}.pack \
  -a -A --progress\
  -t 20 -z -c 50 -C 10000000\
  -s ${id} \
  > vg_vcf/${id}.vcf

bgzip vg_vcf/${id}.vcf
bcftools index vg_vcf/${id}.vcf.gz

