#!/bin/bash 
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=vg-to-calls
#SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/logs/vg_to_calls-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/logs/vg_to_calls-%j.out

module load cactus/3.1.4
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus

export TMPDIR="/scratch/general/nfs1/u6071015/vg_tmp"

PAIR="cactusStripe_TcrGSH2_TcrGUSH2_DTv2"
GENOME1="t_crist_hwy154_cen4119_hap2.fasta.masked"
GENOME2="t_crist_hwy154_cen4280_hap2.fasta.masked"
SAMP1="TcrGSH2"
SAMP2="TcrGUSH2"

#hal2vg ${PAIR}.hal --hdf5InMemory --chop 32 --progress > ${PAIR}.vg

#vg index ${PAIR}.vg -x ${PAIR}.xg -L

#vg snarls ${PAIR}.xg > ${PAIR}.snarls

#vg convert -f ${PAIR}.vg > ${PAIR}.gfa

#vg gbwt --num-jobs 16 --gbz-format -g ${PAIR}.gbz -G ${PAIR}.gfa

vg prune -r -p -t 16 ${PAIR}.vg > ${PAIR}.pruned.vg

vg index ${PAIR}.pruned.vg \
         -L -j ${PAIR}.pruned.dist 

#vg giraffe -b hifi -Z ${PAIR}.gbz \
#	-f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/${GENOME1} \
#	-p > ${SAMP1}.gam

#vg giraffe -b hifi -Z ${PAIR}.gbz \
#	-f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/${GENOME2} \
#	-p > ${SAMP2}.gam

#vg pack ${PAIR}.vg \
#        -g ${SAMP1}.gam \
#        -g ${SAMP2}.gam \
#        -o ${PAIR}.pack

#vg call -A -c 50 -r ${PAIR}.snarls \
#	--threads 6 -S ${SAMP1} \
#	-k ${PAIR}.pack \
#	${PAIR}.vg > ${PAIR}.vcf.gz
