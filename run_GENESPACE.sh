#!/bin/bash 
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos=gompert-grn
#SBATCH --job-name=GENESPACE
#SBATCH --error=/scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/GENESPACE-%j.err
#SBATCH --output=/scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/GENESPACE-%j.out

#load modules
module load orthofinder
module load R

cd /scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/

echo "start GENESPACE"

Rscript genespace_TIMEMA.R

echo "GENESPACE done"
