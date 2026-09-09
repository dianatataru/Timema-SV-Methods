#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/glgenest_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/glgenest_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=glgenest
#SBATCH --qos gompert-grn

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf

## posterior mode
./gl2genestMax af_REF_all.txt  REF_all.gl
## posterior mean
./gl2genest af_REF_all.txt   REF_all.gl
