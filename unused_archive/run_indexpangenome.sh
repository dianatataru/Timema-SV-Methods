#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/indexpangenome_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/indexpangenome_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=indexpangenome
#SBATCH --qos gompert-grn
#SBATCH --mem=200G

### LOAD MODULES ###
module load bwa

### ASSIGN VARIABLES  ###
pangenome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.sv.gfa.fa.gz"
echo "pangenome=$pangenome"

#index fasta
bwa index ${pangenome}
