#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=cactus
##SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/logs/cactus-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/logs/cactus-%j.out

#cd /scratch/general/nfs1/u6071015/cactusNp

module load cactus/2.7.2

#cactus timemajobStore_TcrGSH2_TcrGSR1v2  /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/config/cactusTcrGSH2_TcrGSR1.txt cactusStripe_TcrGSH2_TcrGSR1v2.hal --maxCores 80
#cactus timemajobStore_TcrGSH2_TcrGUSH2v2 --restart /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/config/cactusTcrGSH2_TcrGUSH2.txt cactusStripe_TcrGSH2_TcrGUSH2_DTv2.hal --maxCores 80
#cactus timemajobStore_TcrGSH2_TcrGUSH1v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/config/cactusTcrGSH2_TcrGUSH1.txt cactusStripe_TcrGSH2_TcrGUSH1_DTv2.hal --maxCores 80
#cactus timemajobStore_TcrGSH2_TcrGSR2v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/config/cactusTcrGSH2_TcrGSR2.txt cactusStripe_TcrGSH2_TcrGSR2_DTv2.hal --maxCores 80
#cactus timemajobStore_TcrGSH2_TcrGSH1v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/config/cactusTcrGSH2_TcrGSH1.txt cactusStripe_TcrGSH2_TcrGSH1_DTv2.hal --maxCores 80
#cactus timemajobStore_TcrGSH2_TcrGUSR1v2 --restart /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/config/cactusTcrGSH2_TcrGUSR1.txt cactusStripe_TcrGSH2_TcrGUSR1_DTv2.hal --maxCores 80
#cactus timemajobStore_TcrGSH2_TcrGUSR2v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/config/cactusTcrGSH2_TcrGUSR2.txt cactusStripe_TcrGSH2_TcrGUSR2_DTv2.hal --maxCores 80

#cp /scratch/general/nfs1/u6071015/cactusNp/cactusStripe_TcrGSH2_TcrGSR1v2.hal /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusStripe_TcrGSH2_TcrGSR1v2.hal
#cp /scratch/general/nfs1/u6071015/cactusNp/cactusStripe_TcrGSH2_TcrGUSH2_DTv2.hal /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusStripe_TcrGSH2_TcrGUSH2_DTv2.hal
#cp /scratch/general/nfs1/u6071015/cactusNp/cactusStripe_TcrGSH2_TcrGUSH1_DTv2.hal /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusStripe_TcrGSH2_TcrGUSH1_DTv2.hal
#cp /scratch/general/nfs1/u6071015/cactusNp/*_DTv2.hal /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus

#Summarize Mutations
halSummarizeMutations cactusStripe_TcrGSH2_TcrGSH1_DTv2.hal
#halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSR2_DTv2.hal
#halSummarizeMutations cactusStripe_TcrGSH2_TcrGSR1_DTv2.hal
#halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSH2_DTv2.hal 
#halSummarizeMutations cactusStripe_TcrGSH2_TcrGSR2_DTv2.hal
#halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSR1_DTv2.hal
#halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSH1_DTv2.hal
