# Timema-SV-Methods
Identifying how different methods of structural variant (SV) detection identify SVs, and whether Does SV size, frequency, or age have an effect on SV detection in different methods? Using data from this Gompert et al. 2025 (https://www.science.org/doi/10.1126/science.adp3745), and specifically focused on inversions and translocations.

## Pangenome Creation
Starting off with making pangenomes with 1) 4 hwy154 genomes and 2) 8 genomes (4 hwy154 and 4 refugio) in cactus minigraph (paper:https://www.nature.com/articles/s41587-023-01793-w, documentation: https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).

|  POP  | STRIPE |   ID  |
|-------|--------|-------|
|HWY154 |striped |cen4119|
|HWY154 |  green |cen4280|
|REFUGIO|striped |cen4122|
|REFUGIO|  green |cen4120|

Working directory can be found in ```/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/```. I copied the soft masked genomes to subdirectory "genomes" in this directory. Note, Minigraph-Cactus ignores softmasking. 

```
ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gs_hap_cen4119/HiRise/Hap1/chroms_final_assembly.fasta.masked t_crist_hwy154_cen4119_hap1.fasta.masked 
ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gs_hap_cen4119/HiRise/Hap2/chroms_final_assembly.fasta.masked t_crist_hwy154_cen4119_hap2.fasta.masked 
ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap1/chroms_final_assembly.fasta.masked t_crist_hwy154_cen4280_hap1.fasta.masked
ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap2/chroms_final_assembly.fasta.masked t_crist_hwy154_cen4280_hap2.fasta.masked
ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/chroms_final_assembly.fasta.masked t_crist_refug_cen4122_hap1.fasta.masked
ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap2/chroms_final_assembly.fasta.masked t_crist_refug_cen4122_hap2.fasta.masked
ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap1/chroms_final_assembly.fasta.masked t_crist_refug_cen4120_hap1.fasta.masked
ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap2/chroms_final_assembly.fasta.masked t_crist_refug_cen4120_hap2.fasta.masked
```
Make the HWY154.txt input file (reference can't start with same name as other samples):

```
Hap1_t_crist_hwy154_cen4119 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap1.fasta.masked
t_crist_hwy154_cen4119.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked
t_crist_hwy154_cen4280.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap1.fasta.masked
t_crist_hwy154_cen4280.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked 
```
Make the HWY154_REF.txt input file:

```
t_crist_hwy154_cen4119.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap1.fasta.masked
t_crist_hwy154_cen4119.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked
t_crist_hwy154_cen4280.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap1.fasta.masked
t_crist_hwy154_cen4280.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked 
t_crist_refug_cen4122.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4122_hap1.fasta.masked
t_crist_refug_cen4122.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4122_hap2.fasta.masked
t_crist_refug_cen4120.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4120_hap1.fasta.masked   
t_crist_refug_cen4120.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4120_hap2.fasta.masked
```

Cactus minigraph does require designating a reference, and you can actually designate multiple to use as the basis for VCF files (with command --vcfReference). Other commands that I haven't run yet but could include if the graph seems weird, are --permissiveContigFilter and --noSplit (which disables chromosome splitting). Now starting the minigraph pipeline with script ```run_cactus-pangenome.sh```:

Not working in SBATCH. error is: environment:ancestorsmlmp.py not found. 

```
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus-pangenome
#SBATCH --qos=gompert-grn
#SBATCH -e /scratch/general/nfs1/u6071015/cactusNp/timema/cactus-pangenome-%j.err
#SBATCH -o /scratch/general/nfs1/u6071015/cactusNp/timema/cactus-pangenome-%j.out

module load cactus/3.0.1

cd /scratch/general/nfs1/u6071015/cactusNp/timema/

cactus-pangenome timemaJS \
  /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/HWY154.txt \
  --outDir /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus \
  --outName HWY154 \
  --reference Hap1_t_crist_hwy154_cen4119 \
  --maxCores 24 \
  --vcf --giraffe --gfa --gbz --viz

```
Also works as an interactive script:

```
salloc --time=10:00:00 --ntasks 24 --nodes=1 --account=gompert --partition=gompert-grn --qos=gompert-grn --mem=100G
cd /scratch/general/nfs1/u6071015/cactusNp/timema/
module load cactus/3.0.1
#module load apptainer/1.4.0  
#APPTAINERENV_PREPEND_PATH="/home/cactus/bin"

cactus-pangenome timemaJS \
  /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/HWY154.txt \
  --outDir /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus \
  --outName HWY154 \
  --reference Hap1_t_crist_hwy154_cen4119 \
  --maxCores 12 \
  --vcf --giraffe --gfa --gbz
```
Even with 100G, past the limit (107G), increase for 8 genomes. Getting this error:

```
Got message from job at time 11-10-2025 17:18:07: Job used more disk than requested. For CWL, consider increasing the outdirMin requirement, otherwise, consider increasing the disk requirement. Job 'unzip_gz' kind-unzip_gz/instance-b62o1k6q v1 used 102.02% disk (1.9 GiB [2029998080B] used, 1.9 GiB [1989727325B] requested).
```
Finished after about 10 hours in interactive job.

Quickly visualize output vg easily with tube map onlione demo:https://vgteam.github.io/sequenceTubeMap/. but the file has to be under 5 mb so subset or run on cluster.
Possible analysis can be done in odgi: https://odgi.readthedocs.io/en/latest/

### Investigating Cactus Output

Output with Hap1_t_crist_hwy154_cen4119 as reference:
```
halStats HWY154.full.hal
INFO:    gocryptfs not found, will not be able to use gocryptfs

hal v2.2
(t_crist_hwy154_cen4280.1:1,t_crist_hwy154_cen4119.2:1,t_crist_hwy154_cen4280.2:1,Hap1_t_crist_hwy154_cen4119:1)Anc0;

GenomeName,         NumChildren, Length,   NumSequences, NumTopSegments, NumBottomSegments
Anc0,                     4,     1697544243, 645100,     0,             14123991
t_crist_hwy154_cen4280.1, 0,     1204896739, 13,         10845222,       0
t_crist_hwy154_cen4119.2, 0,     1226560494, 13,         10915768,       0
t_crist_hwy154_cen4280.2, 0,     1215314917, 13,         10916299,       0
Hap1_t_crist_hwy154_cen4119, 0,  1220429573, 13,         10949670,       0
```

how the raw.vcf file looks:
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	t_crist_hwy154_cen4119	t_crist_hwy154_cen4280
Scaffold_10__1_contigs__length_74320458	8053	>36>38	A	ATA	60	AC=1;AF=1;AN=1;AT=>36>38,>36>37>38;NS=1;LV=0	GT	.|.	1|.
Scaffold_10__1_contigs__length_74320458	13799	>43>45	T	TC	60	AC=1;AF=1;AN=1;AT=>43>45,>43>44>45;NS=1;LV=0	GT	.|.	1|.
Scaffold_10__1_contigs__length_74320458	14013	>45>47	GA	G	60	AC=1;AF=1;AN=1;AT=>45>46>47,>45>47;NS=1;LV=0	GT	.|.	1|.
```

