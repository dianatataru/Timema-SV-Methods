#!/bin/bash 
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus-syn
#SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/cactus-syn-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/cactus-syn-%j.out

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/

module load cactus/2.7.2

#halSynteny --queryGenome TcrGSR2 --targetGenome TcrGSH2 cactusStripe_TcrGSH2_TcrGSR2_DTv2.hal cactusStripe_TcrGSH2_TcrGSR2.psl
halSynteny --queryGenome TcrGSR1 --targetGenome TcrGSH2 cactusStripe_TcrGSH2_TcrGSR1_DTv2.hal cactusStripe_TcrGSH2_TcrGSR1.psl
halSynteny --queryGenome TcrGSH1 --targetGenome TcrGSH2 cactusStripe_TcrGSH2_TcrGSH1_DTv2.hal cactusStripe_TcrGSH2_TcrGSH1.psl
halSynteny --queryGenome TcrGUSH1 --targetGenome TcrGSH2 cactusStripe_TcrGSH2_TcrGUSH1_DTv2.hal cactusStripe_TcrGSH2_TcrGUSH1.psl
halSynteny --queryGenome TcrGUSH2 --targetGenome TcrGSH2 cactusStripe_TcrGSH2_TcrGUSH2_DTv2.hal cactusStripe_TcrGSH2_TcrGUSH2.psl
halSynteny --queryGenome TcrGUSR2 --targetGenome TcrGSH2 cactusStripe_TcrGSH2_TcrGUSR2_DTv2.hal cactusStripe_TcrGSH2_TcrGUSR2.psl
halSynteny --queryGenome TcrGUSR1 --targetGenome TcrGSH2 cactusStripe_TcrGSH2_TcrGUSR1_DTv2.hal cactusStripe_TcrGSH2_TcrGUSR1.psl

