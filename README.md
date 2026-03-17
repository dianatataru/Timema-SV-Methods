# Timema-SV-Methods
Identifying how different methods of structural variant (SV) detection identify SVs, and whether Does SV size, frequency, or age have an effect on SV detection in different methods? Using data from this Gompert et al. 2025 (https://www.science.org/doi/10.1126/science.adp3745), and specifically focused on inversions and translocations. Here is some helpful background on pangenomics: https://pangenome.github.io/

## Pangenome Creation
Starting off with making pangenomes with 1) 4 hwy154 genomes and 2) 8 genomes (4 hwy154 and 4 refugio) in cactus minigraph (paper:https://www.nature.com/articles/s41587-023-01793-w, documentation: https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).

|  POP  | STRIPE |   ID  | SHORT |HAP|    SCIENCE BOUNDS    |
|-------|--------|-------|-------|---|----------------------|
|HWY154 |striped |cen4119| H GS  | 2 |24,457,103–39,030,359 |
|HWY154 |  green |cen4280| H GUS | 1 |24,803,527–44,121,870 |
|REFUGIO|striped |cen4122| R GS  | 1 |22,442,098–65,729,704 |
|REFUGIO|  green |cen4120| R GUS | 1 |22,220,178–65,829,835 |

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
t_crist_hwy154_cen4119.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap1.fasta.masked
Hap2_t_crist_hwy154_cen4119.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked
t_crist_hwy154_cen4280.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap1.fasta.masked
t_crist_hwy154_cen4280.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked 
```
Make the HWY154_REF.txt input file:

```
t_crist_hwy154_cen4119.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap1.fasta.masked
Hap2_t_crist_hwy154_cen4119.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked
t_crist_hwy154_cen4280.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap1.fasta.masked
t_crist_hwy154_cen4280.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked 
t_crist_refug_cen4122.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4122_hap1.fasta.masked
t_crist_refug_cen4122.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4122_hap2.fasta.masked
t_crist_refug_cen4120.1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4120_hap1.fasta.masked   
t_crist_refug_cen4120.2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4120_hap2.fasta.masked
```

Cactus minigraph does require designating a reference, and you can actually designate multiple to use as the basis for VCF files (with command --vcfReference). The cactus publication runs on twp different references and chooses the longest one. Other commands that I haven't run yet but could include if the graph seems weird, are --permissiveContigFilter and --noSplit (which disables chromosome splitting). Also, --vcfbub by default flattens the vcf and removes nested variants. If I want that not to happen, I have to set --vcfbub 0. This might make downstream annotation harder. Other things to note on this graph: minigraph only uses SVs > 50 bp in graph construction, and also clips out stretches of sequences >= 10 kb that do not align to minigraph. Now starting the minigraph pipeline with script ```run_cactus-pangenome.sh```:

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
  --vcfbub 0 --giraffe --gfa --gbz --viz --chrom-og

cactus-pangenome timema8hapJS \
  /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/HWY154_REF.txt \
  --outDir /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus \
  --outName HWY154_4119Hap2 \
  --reference Hap2_t_crist_hwy154_cen4119 \
  --maxCores 24 \
  --vcfbub 0 --chrom-og --viz

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

Quickly visualize output vg with tube map onlione demo:https://vgteam.github.io/sequenceTubeMap/. but the file has to be under 5 mb so subset or run on local computer.


### Investigating Cactus Pangenome Output

halStats output:

```
#with Hap1_t_crist_hwy154_cen4119 as reference
GenomeName,         NumChildren, Length,   NumSequences, NumTopSegments, NumBottomSegments
Anc0,                     4,     1697544243, 645100,     0,             14123991
t_crist_hwy154_cen4280.1, 0,     1204896739, 13,         10845222,       0
t_crist_hwy154_cen4119.2, 0,     1226560494, 13,         10915768,       0
t_crist_hwy154_cen4280.2, 0,     1215314917, 13,         10916299,       0
Hap1_t_crist_hwy154_cen4119, 0,  1220429573, 13,         10949670,       0

#with Hap2_t_crist_hwy154_cen4119 as reference
GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
Anc0, 4, 1682922426, 645056, 0, 14164284
t_crist_hwy154_cen4280.2, 0, 1215314917, 13, 10944330, 0
t_crist_hwy154_cen4280.1, 0, 1204896739, 13, 10883086, 0
Hap2_t_crist_hwy154_cen4119.2, 0, 1226560494, 13, 11005277, 0
t_crist_hwy154_cen4119.1, 0, 1220429573, 13, 10903134, 0

```
and for pangenomes with Refugio included:
```
#Hwy 154 Striped Haplotype 1 as Reference: 
GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
Anc0, 8, 2317835600, 1465030, 0, 24034530
t_crist_refug_cen4120.2, 0, 913092304, 11, 13453958, 0
t_crist_refug_cen4122.2, 0, 918956906, 11, 13528212, 0
t_crist_refug_cen4120.1, 0, 919573812, 11, 13420062, 0
t_crist_refug_cen4122.1, 0, 1227621598, 13, 17085224, 0
t_crist_hwy154_cen4280.1, 0, 1204896739, 13, 17105678, 0
t_crist_hwy154_cen4119.2, 0, 1226560494, 13, 17235926, 0
t_crist_hwy154_cen4280.2, 0, 1215314917, 13, 17244420, 0
Hap1_t_crist_hwy154_cen4119, 0, 1220429573, 13, 17266391, 0

#Hwy 154 Striped Haplotype 2 as Reference:
GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
Anc0, 8, 2597813988, 1465202, 0, 24632873
t_crist_refug_cen4120.1, 0, 1239970858, 12, 15785477, 0
t_crist_refug_cen4122.1, 0, 1227621598, 13, 17528880, 0
t_crist_hwy154_cen4280.2, 0, 1215314917, 13, 17664031, 0
t_crist_hwy154_cen4280.1, 0, 1204896739, 13, 17523127, 0
t_crist_refug_cen4122.2, 0, 918956906, 11, 13595387, 0
Hap2_t_crist_hwy154_cen4119.2, 0, 1226560494, 13, 17722387, 0
t_crist_refug_cen4120.2, 0, 1235469353, 12, 15808690, 0
t_crist_hwy154_cen4119.1, 0, 1220429573, 13, 17566620, 0

```
The outputs of HalSummarizeMutations are in this google sheet:https://docs.google.com/spreadsheets/d/1sTRpJKJHh38i-38SDlRJKjfCCZoWqW8oGsbMvjLeViY/edit?usp=sharing
This is how the raw.vcf file looks (head):
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	t_crist_hwy154_cen4119	t_crist_hwy154_cen4280
Scaffold_10__1_contigs__length_74320458	8053	>36>38	A	ATA	60	AC=1;AF=1;AN=1;AT=>36>38,>36>37>38;NS=1;LV=0	GT	.|.	1|.
Scaffold_10__1_contigs__length_74320458	13799	>43>45	T	TC	60	AC=1;AF=1;AN=1;AT=>43>45,>43>44>45;NS=1;LV=0	GT	.|.	1|.
Scaffold_10__1_contigs__length_74320458	14013	>45>47	GA	G	60	AC=1;AF=1;AN=1;AT=>45>46>47,>45>47;NS=1;LV=0	GT	.|.	1|.
```
From Science paper (Gompert et al. 2025), length of chromosomes:

Table S1: Homologous relationships among chromosome-size scaffolds for the T. cristinae
genomes. Chromosome 13 is the X sex chromosome. Abbreviations are as follows: Chr = chromosome,
U GS = unphased stripe genome from (82), R = Refugio, H = Hwy154, GS1 = striped
haplotype 1, GS2 = striped haplotype 2, GUS1 = green haplotype 1, and GUS2 = green haplotype.
4119, 4122=striped and 4120,4280=green. *reference for these analyses is H GS2*

|Chr |U GS |R GS1 |R GS2 |R GUS1 |R GUS2 |H GS1 |H GS2 |H GUS1 |H GUS2|
|----|-----|------|------|-------|------ |------|------|-------|------|
|  1 | 8483|  12  |  10  |   4   |   6   |  13  |  13  |  22   |  15  |
| 2  |14640|   6  |   8  |  11   |   9   |   5  |   6  |  23   |   1  |
| 3  |42935|   2  |   1  |   1   |   1   |   3  |   2  |  16   |   3  |
| 4  |42912|   1  |   1  |   1   |   1   |   1  |   1  |  64   |  35  |
| 5  |18722|   7  |   2  |   6   |   4   |  12  |  12  |   5   |  10  |
| 6  |9928 |   8  |   5  |   7   |   5   |   4  |   5  |  11   |  44  |
| 7  |10660|  10  |   7  |  10   |  11   |  10  |   8  |  54   |   7  |
| 8  |7748 |  11  |   9  |   3   |   3   |  11  |   4  |   7   |  23  |
| 9  |16151|   5  |   4  |   9   |   8   |   8  |   9  |  46   |  21  |
| 10 |14160|   4  |   3  |   8   |  10   |   7  |   7  |  15   |  16  |
| 11 |12033|   9  |   6  |   5   |   7   |   9  |  10  |   2   |  12  |
| 12 |12380|  13  |  12  |  12   |  12   |   6  |  11  |   1   |  36  |
| 13 |14101|   3  |  11  |   2   |   2   |   2  |   3  |  36   |   8  |

### Using sequenceTubeMap

Downloaded to local computer using these instructions: https://github.com/vgteam/sequenceTubeMap?tab=readme-ov-file

Then to run tube map on local computer in terminal:

```
cd ~/Desktop/GitHub/sequenceTubeMap
nvm use
npm run serve
```
I uploaded the .gbz and .gaf file from the cluster into the folder ~/Desktop/GitHub/sequenceTubeMap/exampleData following these instructions: https://github.com/vgteam/sequenceTubeMap/blob/master/doc/data.md. I went here to visualize: http://localhost:3000.

The .gbz gets mounted as the graph and haplotype, while the .gaf can be mounted as the reads. Make sure to index the .gaf file with tabix (htslib), and upload the index file to /exampleData/ as well. 

### Using ODGI

Downloaded using conda on kingspeak, because granite seems to not be working at this moment.

```
module load miniforge3
conda create -n odgi
conda activate odgi
mamba install -c bioconda odgi
mamba install -c bioconda vg

```

need to convert .vg to .gfa file to .og format. Sorting in OG format will also help with the complexity found in the current .viz graphs. If we want to make the loopy line plots, my understanding is that you have to do the following. Note, you have to run this chromosome by chromsome. I went into the chrom_alignments folder and moved older alignments into old_alignments subdirectory to run all of this on the HWY154_4119Hap2 reference with all 8 genomes. 

Messing around with cactus outoput files:
```
#start interactive job
salloc --time=06:00:00 --ntasks 24 --nodes=1 --account=gompert-kp --partition=gompert-kp --mem=100G

module load miniforge3
conda activate odgi

#or confert .vg to v1.0 gfa
vg convert -f Scaffold_9__2_contigs__length_79556474.vg  -W > Scaffold_9__2_contigs__length_79556474.gfa

#creat .og file
odgi build -g  Scaffold_9__2_contigs__length_79556474.gfa  -o Scaffold_9__2_contigs__length_79556474.og

#sort .og file (-Y selects the PG-SGD algorithm for sorting, many options to tweak this)
# I didnt change the max number of iterations here (default 30) but I could using -x
odgi sort -i Scaffold_9__2_contigs__length_79556474.og --threads 20 -P -C /scratch/general/nfs1/u6071015/odgi -o Scaffold_9__2_contigs__length_79556474_sorted.og

#visualize sorted graph (this creates a png that is too big to view)
odgi viz -i Scaffold_9__2_contigs__length_79556474_sorted.og -o Scaffold_9__2_contigs__length_79556474_sorted.svg -x 5000 -y 2000

#SV calling with vg
vg snarls Scaffold_10__2_contigs__length_75648701.vg > Scaffold_10__2_contigs__length_75648701.snarls

#this hasn't resulted in any sort of interpretable output, now following some of the analysis in here: https://cpang.netlify.app/post/day-3-bacterial-pangenomics/
#index in vg
vg index -x Scaffold_9__2_contigs__length_79556474.xg Scaffold_9__2_contigs__length_79556474.vg

#visualize in vg
vg viz -x Scaffold_9__2_contigs__length_79556474.xg -o Scaffold_9__2_contigs__length_79556474.svg

```
Making the loopy odgi draw graphs with odgi_draw.sh:
```
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=odgi-draw
#SBATCH -e odgi-draw-%j.err
#SBATCH -o odgi-draw-%j.out

module load cactus/3.0.1

cactus-graphmap-join /scratch/general/nfs1/u6071015/cactusNp/timema/timema8hapJS \
   --vg /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus/chrom-alignments/*.vg \
   --reference Hap2_t_crist_hwy154_cen4119 \
   --outDir /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus/HWY154_REF_4119Hap2/chroms \
   --outName HWY154_REF_4119Hap2 --draw

```
making the sorted viz graphs (hopefully less paths, more interpretable) with run_odgi_viz.sh:

```
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=odgi-viz
#SBATCH -e odgi-viz-%j.err
#SBATCH -o odgi-viz-%j.out

module load miniforge3
conda activate odgi

WORKDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments"
cd $WORKDIR

SCAFFOLDS=(
  Scaffold_1__1_contigs__length_160647932
  Scaffold_2__1_contigs__length_157594471
  Scaffold_3__2_contigs__length_137956696
  Scaffold_4__1_contigs__length_97222829
  Scaffold_5__1_contigs__length_83128659
  Scaffold_6__1_contigs__length_78844258
  Scaffold_7__1_contigs__length_75018798
  Scaffold_8__1_contigs__length_71271319
  Scaffold_9__2_contigs__length_79556474
  Scaffold_10__2_contigs__length_75648701
  Scaffold_11__2_contigs__length_80009992
  Scaffold_13__3_contigs__length_82050896
)

for SCAFFOLD in "${SCAFFOLDS[@]}"; do
    echo "Processing ${SCAFFOLD}"


odgi build -g ${SCAFFOLD}.gfa  -o ${SCAFFOLD}.og -O

odgi sort -i ${SCAFFOLD}.og --threads 20 -P -Y -C /scratch/general/nfs1/u6071015/odgi -o ${SCAFFOLD}_sortedPGSGD.og

#make a file that lists just the genome paths. the graph will be too big with all MINIGRAPH paths.
mapfile -t TCRIST_PATHS < <(
  odgi paths -L -i ${SCAFFOLD}_sortedPGSGD.og | grep 't_crist'
)
PATHFILE="display_paths_${SCAFFOLD}.txt"
printf "%s\n" "${TCRIST_PATHS[@]}" > ${PATHFILE}

odgi viz \
    -i ${SCAFFOLD}_sortedPGSGD.og \
    -o ${SCAFFOLD}_sortedPGSGD.png \
	-p ${PATHFILE}
done
```


### Using Pantree to describe SVs
Following this preprint: https://www.biorxiv.org/content/10.1101/2025.08.04.668502v1

```
pwd /uufs/chpc.utah.edu/common/home/u6071015/software
module load miniforge3
pip install uv
git clone ssh:://git@github.com/oclb/graph_var.git
cd graph_var
uv venv
#Using CPython 3.13.9
#Creating virtual environment at: .venv
#Activate with: source .venv/bin/activate
uv sync
# Built pantree @ file:///uufs/chpc.utah.edu/common/home/u6071015/software/pantree
```
Okay, so now I am going to move to the /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree folder
and here is the sbatch script to run the python script, which I will call run_pantree.sh:
```
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=pantree
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/pantree-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/pantree-%j.out

SCAFF="Scaffold_2__1_contigs__length_157594471"

module load miniforge3

#convert vg
conda activate odgi
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments
vg convert -f  ${SCAFF}.vg -W > ${SCAFF}.gfa
conda deactivate

#run pantree
cd /uufs/chpc.utah.edu/common/home/u6071015/software/pantree
source .venv/bin/activate
export PYTHONPATH="/uufs/chpc.utah.edu/common/home/u6071015/software/pantree:${PYTHONPATH}"

uv run pantree /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/${SCAFF}.gfa \
        /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/${SCAFF}_pantree.vcf.gz \
        --ref-name Hap2_t_crist_hwy154_cen4119.2 \
        --log-path /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/${SCAFF}.log \
        --chr-id ${SCAFF} \
        --priority-samples t_crist_hwy154_cen4119.1,t_crist_hwy154_cen4280.1,t_crist_hwy154_cen4280.2,t_crist_refug_cen4122.1,t_crist_refug_c$

```
My old run timed out due to some bugs in the program. I changed the name of the old program to pantree_OLD and downloaded the updated program in /uufs/chpc.utah.edu/common/home/u6071015/software/pantree. It had an OOM killed event trying to run the entire genome, so I have to run it scaffold by scaffold.

The output header of the pantree.vcf.gz looks like this:

```
##fileformat=VCFv4.2
##source=pantree v0.2.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 1 if ALT is present, 0 if absent, . if missing">
##FORMAT=<ID=CR,Number=R,Type=Integer,Description="Number of times visiting the REF allele">
##FORMAT=<ID=CA,Number=A,Type=Integer,Description="Number of times visiting the ALT allele">
##INFO=<ID=NR,Number=1,Type=String,Description="Non-reference allele">
##INFO=<ID=VT,Number=1,Type=String,Description="Variant type">
##INFO=<ID=TP,Number=1,Type=Integer,Description="Tree position of the variant edge's branch point">
##INFO=<ID=RC,Number=1,Type=Integer,Description="The REF allele count">
##INFO=<ID=AC,Number=A,Type=Integer,Description="The ALT allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=HP,Number=.,Type=String,Description="Haplotype positions at reference tree edge (haplotype:position)">
##INFO=<ID=TR_MOTIF,Number=1,Type=String,Description="Tandem repeat motif">
##INFO=<ID=NIA,Number=1,Type=Integer,Description="Nearly identical alleles (1=yes, 0=no)">
##INFO=<ID=UIDX,Number=1,Type=Integer,Description="Index of node u">
```
Now to summarize the output vcf using code from their paper (https://github.com/ShenghanZhang1123/graph_var_analysis/blob/main/notebooks/generating_data_analysis.ipynb) in a script I wrote called pantree_summary.py:

Run it using this sbatch script, it outputs into the ```summary``` directory:

```
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
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
SCAFF="Scaffold_9__2_contigs__length_79556474"

#mkdir summary

#summarize SVs (edited from pantree manuscript, puts output in /summary subdir of working directory)
python pantree_summary.py --vcf /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/${SCAFF}_pantree.vcf.gz --chrom ${SCAFF}

```
The output from this summary has only non-linear variants. I'm wondering if I also have to include a sample list for the --ref-name to be used? I'm going to try this again on Scaffold 12, the shortest one, with the samples listed. I changed the name of the original pantree vcf made to Scaffold_12__1_contigs__length_47609450_nonlinear_pantree.vcf.gz, so that it is not overwritten. This is what I added to the end of the pantree run line: 
--priority-samples t_crist_hwy154_cen4119.1,t_crist_hwy154_cen4280.1,t_crist_hwy154_cen4280.2,t_crist_refug_cen4122.1,t_crist_refug_cen4122.2,t_crist_refug_cen4120.1,t_crist_refug_cen4120.2

I had to change the run time to 48 hours due to some maintenance on the cluster, should turn it back to 7 days when I can to run the remain scaffolds (2-8, 10, 11)

To just output inversion vcf:

```
salloc --time=06:00:00 --ntasks 12 --nodes=1 --account=gompert --partition=gompert-grn --qos gompert-grn
zcat Scaffold_9__2_contigs__length_79556474_pantree.vcf.gz \
| awk '
  /^#/ { print; next }
  $8 ~ /(^|;)VT=INV(;|$)/
' \
| gzip > Scaffold_9__2_contigs__length_79556474_pantree_inversions_only.vcf.gz
```

For Scaffold 9, there are 3 inversions (found by my old summary script and then also by running   
```
 bcftools query -f '%POS\n' Scaffold_9__2_contigs__length_79556474_pantree_inversions_only.vcf.gz
 ```
Wierdly, it says that the positions for all of these inversions are 1, which doesn't make sense.
Also, wierdly, one of the inversions has no Ref or alt? The second inversion is huge.

to subset to SVs no SNPs or bigger than 50 bp:

```
salloc --time=06:00:00 --ntasks 12 --nodes=1 --account=gompert --partition=gompert-grn --qos gompert-grn
module load bcftools
bcftools query \
  -f '%CHROM\t%POS\t%ID\t%INFO/AC\t%INFO/RC\t%INFO/VT\t%INFO/TP\n' \
  Scaffold_4__1_contigs__length_97222829_pantree.vcf.gz \
  > Scaffold_4__1_contigs__length_97222829_pantree.vcf_minimal.tsv

 awk -F'\t' '$6 != "SNP"' Scaffold_4__1_contigs__length_97222829_pantree.vcf_minimal.tsv \
  > Scaffold_4_noSNPs.tsv

bcftools view -i 'strlen(INFO/NR) > 50'  Scaffold_4__1_contigs__length_97222829_pantree.vcf.gz -Oz -o Scaffold_4__1_contigs__length_97222829_pantree.len50bpplus.vcf.gz


bcftools query -f '%ID\t%REF\t%ALT\n' Scaffold_4__1_contigs__length_97222829_pantree_inversions_only.vcf.gz > Scaffold_4_inv_alleles.tsv
awk '{print ">"$1"\n"$3}' Scaffold_4_inv_alleles.tsv > Scaffold_4_inv_alt.fa

module load cactus/3.0.1
vg giraffe -Z hs37d5-pangenome.giraffe.gbz -m hs37d5-pangenome.shortread.withzip.min -z hs37d5-pangenome.shortread.zipcodes -d hs37d5-pangenome.dist -f sim.fq >mapped.gam

vg giraffe \
  -Z HWY154_REF_4119Hap2.d2.gbz \
  -d HWY154_REF_4119Hap2.d2.dist \
  -m HWY154_REF_4119Hap2.d2.min \
  -f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/Scaffold_4_inv_alt.fa \
  > Scaffold_4_inv_alt.gam

KILLED

vg surject \
  -x HWY154_REF_4119Hap2.d2.gbz \
  -b Scaffold_4_inv_alt.gam \
  -p Hap2_t_crist_hwy154_cen4119 \
  > Scaffold_4_inv_alt.bam

#trying to just align pantree swith to one genome
module load minimap2


bcftools query -f '%ID\t%REF\t%ALT\n' Scaffold_4__1_contigs__length_97222829_pantree_inversions_only.vcf.gz | \
  awk '$3 != "." {
    id = $1
    gsub(/>/, "fw", id)
    gsub(/</, "rv", id)
    print ">" id "\n" $3
  }' > Scaffold_4_inv_alt_clean.fa

#t_crist_hwy154_cen4119_hap2
minimap2 -cx asm5 \
  /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked \
  Scaffold_4_inv_alt_clean.fa \
  > t_crist_hwy154_cen4119_hap2_scaff4_inversions.paf

sort -k6,6 -k8,8n t_crist_hwy154_cen4119_hap2_scaff4_inversions.paf >t_crist_hwy154_cen4119_hap2_scaff4_inversions.srt.paf            
paftools.js call t_crist_hwy154_cen4119_hap2_scaff4_inversions.srt.paf  > t_crist_hwy154_cen4119_hap2_scaff4_inversions.srt.paf.txt

#t_crist_refug_cen4120_hap2
minimap2 -cx asm5 \
  /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4120_hap2.fasta.masked \
  Scaffold_4_inv_alt_clean.fa \
  >  t_crist_refug_cen4120_hap2_scaff4_inversions.paf

sort -k6,6 -k8,8n t_crist_refug_cen4120_hap2_scaff4_inversions.paf > t_crist_refug_cen4120_hap2_scaff4_inversions.srt.paf            
paftools.js call t_crist_refug_cen4120_hap2_scaff4_inversions.srt.paf   > t_crist_refug_cen4120_hap2_scaff4_inversions.srt.paf.txt

#t_crist_refug_cen4122_hap1.fasta.masked
minimap2 -cx asm5 \
  /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4122_hap1.fasta.masked \
  Scaffold_4_inv_alt_clean.fa \
  >  t_crist_refug_cen4122_hap1_scaff4_inversions.paf

sort -k6,6 -k8,8n t_crist_refug_cen4122_hap1_scaff4_inversions.paf > t_crist_refug_cen4122_hap1_scaff4_inversions.srt.paf            
paftools.js call t_crist_refug_cen4122_hap1_scaff4_inversions.srt.paf > t_crist_refug_cen4122_hap1_scaff4_inversions.srt.paf.txt

```

### Genome Annotation and GENESPACE visualization

We can use the genespace visualization to validate the inversions and translocations found.
Copying over the braker3 annotations from the Science Paper:

```
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation

cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h1/braker/braker.aa t_crist_refug_green_h1.fa #REDOWNLOAD THIS
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_stripe_h1/brakerV1/braker.aa t_crist_refug_stripe_h1.fa #not found, changed dir name
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h2/braker/braker.aa t_crist_h154_green_h2.fa #not found, reran braker but OLD
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h1/braker/braker.aa t_crist_h154_stripe_h1.fa
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h1/brakerV1/braker.aa t_crist_h154_green_h1.fa #not found, changed dir name
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h2/braker/braker.aa t_crist_refug_green_h2.fa
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_stripe_h2/braker/braker.aa t_crist_refug_stripe_h2.fa
cp ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h2/brakerV1/braker.aa t_crist_h154_stripe_h2.fa #not found, changed dir name


## fix format
perl -p -i -e 's/\.t[0-9]//' *fa

grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h1/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_refug_green_h1.bed
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_green_h2/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_refug_green_h2.bed
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_stripe_h1/brakerV1/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_refug_stripe_h1.bed 
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h2/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_h154_green_h2.bed #not made
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h1/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_h154_stripe_h1.bed 
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h1/brakerV1/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_h154_green_h1.bed 
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_refug_stripe_h2/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_refug_stripe_h2.bed 
grep "gene" ~/../gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_stripe_h2/brakerV1/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_h154_stripe_h2.bed 
```

One of the annotations is missing so need to run braker on it

```
#!/bin/bash 
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=braker
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/braker-%j.err
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/braker-%j.out

source ~/.bashrc

ml braker/3.0.8
ml busco

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/t_crist_hyw154_green_h2

## run braker

braker.pl --genome=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked \
	--prot_seq=/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/proteins.fasta \
	--rnaseq_sets_ids=clean_tcr135.17_0003_R,clean_tcr137.17_0006_R,clean_tcr139.17_0012_R,clean_tcr140.17_0015_R,clean_tcr141.17_0019_R,clean_tcr142.17_0043_R,clean_tcr143.17_0045_R,clean_tcr144.17_0049_R,clean_tcr145.17_0051_R,clean_tcr146.17_0057_R,clean_tcr148.17_0062_R,clean_tcr149.17_0065_R,clean_tcr150.17_0067_R,clean_tcr151.17_0070_R,clean_tcr152.17_0074_R,clean_tcr173.17_0075_R,clean_tcr174.17_0081_R,clean_tcr175.17_0082_R \
	--rnaseq_sets_dirs=/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/rna_seq_for_annotations \
	--AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts \
	--AUGUSTUS_CONFIG_PATH=/uufs/chpc.utah.edu/common/home/u6071015/augustus/config \
	--threads=48 --gff3

## run busco, genome and aa
cd braker
## genome
#busco -i /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked -m geno -o busco_genome_out -l insecta_odb10

## amino acids
busco -i braker.aa -m prot -o busco_aa_out -l insecta_odb10

## Augustus amino acids
cd Augustus #had to add this for it to find the input files
busco -i augustus.hints.aa -m prot -o busco_augustus_aa_out -l insecta_odb10
```

then prepare it like the rest of the files:
```
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation
cp /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/t_crist_hyw154_green_h2/braker/braker.aa t_crist_h154_green_h2.fa
perl -p -i -e 's/\.t[0-9]//' t_crist_h154_green_h2.fa
grep "gene" /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/t_crist_h154_green_h2/braker/braker.gff3 | cut -f 1,4,5,9 | perl -p -i -e 's/ID=//' | perl -p -i -e 's/;//' > t_crist_hyw154_green_h2.bed

#there are peptides missing from the bed genes, which means I need to filter the .fa files to match
module load seqkit/2.8.2
SAMPLE="t_crist_h154_green_h2"
cut -f4 ${SAMPLE}.bed | sort -u > bed.ids
seqkit grep -f bed.ids ${SAMPLE}.fa -o ${SAMPLE}_filtered.fa
cp ${SAMPLE}_filtered.fa /scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/peptide/${SAMPLE}.fa
#for samples where the filtered matched unfiltered, I deleted the filtered file.
```

then we can run GENESPACE:

first make folders and move files to those folders:

```
cd /scratch/general/nfs1/u6071015
mkdir GENESPACE_TIMEMA
cd GENESPACE_TIMEMA/
mkdir peptide
mkdir bed
cp /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/*filtered.fa peptide/
cp /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/annotation/*bed bed/

```
then run the program
```
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
```

and this is genespace_TIMEMA.R:

```
#create personal library to write packages to
#dir.create("~/R/x86_64-pc-linux-gnu-library/4.4", recursive = TRUE, showWarnings = FALSE)
#set library paths
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")

#only run devtools download for the first time
#if (!requireNamespace("devtools", quietly = TRUE))
#    install.packages("devtools")
#devtools::install_github("jtlovell/GENESPACE")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("Biostrings", "rtracklayer"))

library(GENESPACE)

#set working directory
wd<-"/scratch/general/nfs1/u6071015/GENESPACE_TIMEMA/"
path2mcscanx<-"~/bin/MCScanX/"

# initalize the run and QC the inputs
gpar<-init_genespace(wd=wd,path2mcscanx=path2mcscanx)

# need to set this
gpar$shellCalls$orthofinder<-"orthofinder"

# accomplish the run
out <- run_genespace(gpar, overwrite = T)

# plot
roi<-data.frame(
		genome=c("t_crist_h154_green_h1","t_crist_h154_green_h2",
			 "t_crist_h154_stripe_h1","t_crist_h154_stripe_h2",
		"t_crist_refug_green_h1","t_crist_refug_green_h2",
    "t_crist_refug_stripe_h1","t_crist_refug_stripe_h2"),
		start=c(0,0,0,0,0,0,0,0,0,0),end=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf))

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))

pdf("syn8way.pdf",width=9,height=5.6)
ripd <- plot_riparian(
	gsParam = out,
	palette = customPal,
        #highlightBed = roi,
	braidAlpha = .3,
	useOrder=TRUE,
	chrFill = "lightgrey",
    addThemes = ggthemes,
	useRegions = FALSE,
	    #invertTheseChrs = invchr,
  	refGenome = "t_crist_h154_stripe_h2",
	genomeIDs = c("t_crist_h154_green_h1","t_crist_h154_green_h2",
			 "t_crist_h154_stripe_h1","t_crist_h154_stripe_h2",
		"t_crist_refug_green_h1","t_crist_refug_green_h2",
    "t_crist_refug_stripe_h1","t_crist_refug_stripe_h2"),
	backgroundColor = NULL)
dev.off()
```


## Pairwise comparison in Progressive Cactus

We are also going to call SVs from the pairwise comparisons, specifically focusing on comparisons between the reference haplotype used for the pangenome (HWY154 Stripe Haplotype2) and the other haplotypes. For this, I am creating softlinks in ''/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus'' to existing hal files, and need to make the hal file for H154 Stripe 2/Refugio Stripe 1 pair. TO do so, I run the script run_cactus.sh in the directory. The script also requires input file cactusTcrGSH2_TcrGSR1.txt:

```
(TcrGSH2:0.010,TcrGSR1:0.010);

TcrGSH2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked
TcrGSR1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4122_hap1.fasta.masked
```
run_cactus.sh:
```
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus
#SBATCH --qos gompert-grn
#SBATCH -e cactus-%j.err
#SBATCH -o cactus-%j.out
#SBATCH --mem=100G

cd /scratch/general/nfs1/u6071015/cactusNp

module load cactus/1.0.0

#missing pair
cactus timemajobStore /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusTcrGSH2_TcrGSR1.txt cactusStripe_TcrGSH2_TcrGSR1.hal --maxCores 80

cp /scratch/general/nfs1/u6071015/cactusNp/cactusStripe_TcrGSH2_TcrGSR1.hal /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusStripe_TcrGSH2_TcrGSR1.hal

#test rerun
cactus timemajobStore_TcrGSH2_TcrGUSH2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusTcrGSH2_TcrGUSH2.txt cactusStripe_TcrGSH2_TcrGUSH2_DT.hal --maxCores 80

cp /scratch/general/nfs1/u6071015/cactusNp/cactusStripe_TcrGSH2_TcrGUSH2_DT.hal /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusStripe_TcrGSH2_TcrGUSH2_DT.hal

module purge
module load cactus/2.7.2
#now running with new cactus version
cactus timemajobStore_TcrGSH2_TcrGSR1v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusTcrGSH2_TcrGSR1.txt cactusStripe_TcrGSH2_TcrGSR1_DTv2.hal --maxCores 80

cactus timemajobStore_TcrGSH2_TcrGSR2v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusTcrGSH2_TcrGSR2.txt cactusStripe_TcrGSH2_TcrGSR2_DTv2.hal --maxCores 80

cactus timemajobStore_TcrGSH2_TcrGSH1v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusTcrGSH2_TcrGSH1.txt cactusStripe_TcrGSH2_TcrGSH1_DTv2.hal --maxCores 80

cactus timemajobStore_TcrGSH2_TcrGUSR1v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusTcrGSH2_TcrGUSR1.txt cactusStripe_TcrGSH2_TcrGUSR1_DTv2.hal --maxCores 80

cactus timemajobStore_TcrGSH2_TcrGUSR2v2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusTcrGSH2_TcrGUSR2.txt cactusStripe_TcrGSH2_TcrGUSR2_DTv2.hal --maxCores 80

cp /scratch/general/nfs1/u6071015/cactusNp/*_DTv2.hal /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus

#Summarize Mutations 
halSummarizeMutations cactusStripe_TcrGSH1_TcrGSH2.hal
halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSH1.hal
halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSR2.hal
halSummarizeMutations cactusStripe_TcrGSH2_TcrGSR1.hal
halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSH2.hal  
halSummarizeMutations cactusStripe_TcrGSH2_TcrGSR2.hal
halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSR1.hal
halSummarizeMutations cactusStripe_TcrGSH2_TcrGUSH2_DT.hal 
```
The output of this ended up with way too many inversions called for my new .hal file (https://docs.google.com/spreadsheets/d/1BqMnLqLyoLgIq9Shhy4W0GRIvbQmTF1AU8gqI19ijA0/edit?gid=0#gid=0), on the scale of 150-250 inversions instead of the normal 5-20 for the existing .hal files. This was because I used Cactus v1 instead of Cactus v2.7.2. I want to use Cactus v2.7.2, because it is better at calling SVs. There is still some slight differences on every cactus run between the same pairs, due to some randomness in the program.

halSummarizeMutations also doesn't come up with identifiers for mutations, so what I really need is a SV caller that develops position/reference-specific (?) identifiers so that I can then tell how many are unique across the pairwise comparisons. The trick is to do this without accidentally just creating another pangenome... it seems like vg might have a way of calling SVs, which would be ideal because I think we are also going to use vg's giraffe-deepvariant workflow for sv calling for GBS data. Also, it seems rigorous compared to other SV callers (Hickey et al. Genome Biology (2020) https://doi.org/10.1186/s13059-020-1941-7).

For the following SV calling, GSH2 is the REF for all. 
```
salloc --time=06:00:00 --ntasks 1 --nodes=1 --account=gompert --partition=gompert-grn --qos gompert-grn --mem=300G
module load cactus/3.0.1
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus

hal2vg cactusStripe_TcrGSH2_TcrGUSH2_DTv2.hal --hdf5InMemory --chop 32 --progress > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg
vg index cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg -L
vg snarls cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls
vg view -j -R cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls.json
wc -l cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls.json
#14234556 cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls.json

vg deconstruct cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg -P TcrGSH2 -e -a > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vcf
bgzip cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vcf
tabix cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vcf.gz

module load bcftools
# subset to SVs greater than 0 bp, throws up error about GT but just rem0ves it
bcftools view -i 'strlen(REF)>50 || strlen(ALT)>50' \
    cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vcf.gz -G -Oz -o cactusStripe_TcrGSH2_TcrGUSH2_min50bp.vcf.gz
tabix -p vcf cactusStripe_TcrGSH2_TcrGUSH2_min50bp.vcf.gz
bcftool stats cactusStripe_TcrGSH2_TcrGUSH2_min50bp.vcf.gz # number of records is 819, I think this is missing a lot of stuff

#truncated due to this VCF parse error:
#Couldn't read GT data: value not a number or '.' at TcrGSH2#0#Scaffold_10__2_contigs__length_75648701:281081 and #TcrGSH2#0#Scaffold_10__2_contigs__length_75648701:2810840

vg chunk \
  -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg \
  -p "TcrGSH2#0#Scaffold_10__2_contigs__length_75648701:2810800-2810900" \
  --snarls cactusStripe_TcrGSH2_TcrGUSH2_DTv2.snarls \
  -g > 2810800-2810900.vg

vg view -d 2810800-2810900.vg > 2810800-2810900.dot
dot -Tpdf 2810800-2810900.dot > 2810800-2810900.pdf

vg paths -v 2810800-2810900.vg -Q TcrGSH2 -L
TcrGSH2#0#Scaffold_10__2_contigs__length_75648701[2810781]
TcrGSH2#0#Scaffold_13__3_contigs__length_82050896[3492648]
TcrGSH2#0#Scaffold_10__2_contigs__length_75648701[6014359]

# this is a translocation, which deconstruct can't handle. I actually want to be using vg call.
# First I need to make associated indexes and a gam for each fasta

vg index -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg \
         -g cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gcsa \
         -L -j cactusStripe_TcrGSH2_TcrGUSH2_DTv2.dist \
         cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg

vg minimizer cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg \
  -d cactusStripe_TcrGSH2_TcrGUSH2_DTv2.dist \
  -o cactusStripe_TcrGSH2_TcrGUSH2_DTv2.min

vg giraffe \
  -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg \
  -g cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gcsa \
  -m cactusStripe_TcrGSH2_TcrGUSH2_DTv2.min \
  -d cactusStripe_TcrGSH2_TcrGUSH2_DTv2.dist \
  -f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked \
  > TcrGSH2.gaf

vg giraffe \
  -x cactusStripe_TcrGSH2_TcrGUSH2_DTv2.xg \
  -g cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gcsa \
  -m cactusStripe_TcrGSH2_TcrGUSH2_DTv2.min \
  -d cactusStripe_TcrGSH2_TcrGUSH2_DTv2.dist \
  -f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4280_hap2.fasta.masked \
  > TcrGUSH2.gaf

vg prune -r -p -t 2 cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg > cactusStripe_TcrGSH2_TcrGUSH2_DTv2.pruned.vg
vg index cactusStripe_TcrGSH2_TcrGUSH2_DTv2.pruned.vg -L -j cactusStripe_TcrGSH2_TcrGUSH2_DTv2.pruned.dist 

vg autoindex --workflow giraffe -g cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gfa \
	-p cactusStripe_TcrGSH2_TcrGUSH2_DTv2 \
	-G cactusStripe_TcrGSH2_TcrGUSH2_DTv2.gbz \
	--threads 2 --target-mem 5G --verbosity 2 -T temp/

# Constructing distance index for Giraffe.
#Killed

vg convert -g TcrGSH2.gaf > TcrGSH2.gam
vg convert -g TcrGUSH2.gaf > TcrGUSH2.gam
```
or as sbatch script 
```
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus
#SBATCH --qos gompert-grn
#SBATCH -e cactus-%j.err
#SBATCH -o cactus-%j.out
#SBATCH --mem=300G

module load cactus/3.0.1
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

vg index ${PAIR}.vg \
         -L -j ${PAIR}.dist 

vg giraffe -b hifi -Z ${PAIR}.gbz \
	-f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/${GENOME1} \
	-p > ${SAMP1}.gam

vg giraffe -b hifi -Z ${PAIR}.gbz \
	-f /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/${GENOME2} \
	-p > ${SAMP2}.gam

vg pack ${PAIR}.vg \
        -g ${SAMP1}.gam \
        -g ${SAMP2}.gam \
        -o ${PAIR}.pack

vg call -A -c 50 -r ${PAIR}.snarls \
	--threads 6 -S ${SAMP1} \
	-k ${PAIR}.pack \
	${PAIR}.vg > ${PAIR}.vcf.gz

```
vg prune worked with a temp dir added, but vg index still oom killing after:

```
[vg prune] Original graph cactusStripe_TcrGSH2_TcrGUSH2_DTv2.vg: 75689771 nodes, 90722944 edges
[vg prune] Built a temporary XG index
[vg prune] Removed all paths
[vg prune] Pruned complex regions: 75689771 nodes, 81196671 edges
[vg prune] Removed small subgraphs: 66856969 nodes, 77460277 edges
Restored graph: 75689771 nodes
[vg prune] Serialized the graph: 75689771 nodes, 90722944 edges
INFO:    gocryptfs not found, will not be able to use gocryptfs
/uufs/chpc.utah.edu/sys/installdir/lmod/8.6-r8/init/bash: line 82: 1265292 Killed                  apptainer exec --nv /uufs/chpc.utah.edu/sys/installdir/cactus/3.1.4/cactus-3.1.4.sif vg $@
slurmstepd: error: Detected 1 oom_kill event in StepId=814509.batch. Some of the step tasks have been OOM Killed.
```

could maybe prune deeper?

```
vg prune \
  -k 16 \
  -X 2 \
  -e 2 \
  -p \
  -t 16 \
  cactusStripe.vg > cactusStripe.pruned.vg
```
New Jay paper does the following with vg deconstruct vcf output:
- ran vcfbub to keep only top-level variant sites (snarls) less than 100 kb in size
- used vcfwave to realign REF and ALT alleles to split nested alleles to separate entries and identify inversions >1kb
- combined vcf files with bcftools concat, added in missing sample coolumns with bcftools query, and used bcftools fixploidy to set allele number for every site and bcftools fill tags to add AF and AC for each each site. Also used bcftools norm to split multiallelic to biallelic

## GBS Data Alignment and Variant Calling from Pangenome with VG

We use the vg-giraffe-pack-call workflow to align and call SVs from the GSH2-8haplotype pangenome(https://link.springer.com/article/10.1186/s13059-020-1941-7#Sec12).

Input graph is the filtered pangenome output form MinigraphCactus that has been vg autoindexed to create the .dist, .shortread.zipcodes, and .shortread.withzip.min indexes. Snarls have been identified using vg snarls. This filter graph is made by removing nodes covered by fewer than 2 haplotypes from the clip graph. The distance index (.dist) is a memory-mapped file. As of vg version 1.48.0, the file will be opened in read+write mode by default. This can cause issues in HPC clusters and other distributed environments, where multiple computers try to access the same distance index file. To avoid this, make the file read-only or use a local copy of the file (chmod 444 HWY154_REF_4119Hap2.d2.dist). Joint variant calling is not possible with vg call. Here, I make a separate vcf for each sample and then will join them all with bcftools merge. I think it is possible to take the output .gam from giraffe and surject it to make bams, and then use GATK gvcf, although it is possible this would return nonsense (https://github.com/vgteam/vg/issues/1416). 

Output from cactus minigraph- generally the clip graph for everything except --giraffe which defaults to the filter graph, and anything odgi-related which defaults to full. 

Make soft linked GBS files:
```
ln -s /uufs/chpc.utah.edu/common/home/gompert-group3/data/sheffield/timema/2013fha_gwas/02_ids_reads/cristinae/2013*bz2 data/
```
now run alignment and variant calling:

```
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=mapGBS
#SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/mapGBS-%A_%a.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/mapGBS-%A_%a.out
#SBATCH --mem=200G
#SBATCH --array=0-1

module load cactus/3.0.1
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS

FILES=(data/*.fq.bz2)
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
id=$(basename "$FILE" .fq.bz2)
echo ID=${id}

PANGENOME_PATH="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2"
PANGENOME="HWY154_REF_4119Hap2.d2"

# Graph alignment (do I need -x graph.xg?
vg giraffe \
  -Z ${PANGENOME_PATH}/${PANGENOME}.gbz \
  -d ${PANGENOME_PATH}/${PANGENOME}.dist \
  -z ${PANGENOME_PATH}/${PANGENOME}.shortread.zipcodes \
  -m ${PANGENOME_PATH}/${PANGENOME}.shortread.withzip.min \
  -f <(bzcat data/${id}.fq.bz2) \
  -t 24 \
  > vg_intermediate/${id}.gam

vg stats -a vg_intermediate/${id}.gam > vg_stats/${id}.stats.txt

# Coverage packing
vg pack \
  -x ${PANGENOME_PATH}/${PANGENOME}.gbz \
  -g vg_intermediate/${id}.gam \
  -Q 5 \
  -t 20 \
  -o vg_intermediate/${id}.pack

# Variant calling
vg call \
  ${PANGENOME_PATH}/${PANGENOME}.gbz \
  -r ${PANGENOME_PATH}/HWY154_REF_4119Hap2.snarls \
  -k vg_intermediate/${id}.pack \
  -a -A --progress\
  -t 20 -z -c 50 -C 10000000\
  -s ${id} \
  > vg_vcf/${id}.vcf

```

output vcf:
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_4__1_contigs__length_97222829,length=97222830>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_9__2_contigs__length_79556474,length=79556476>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_1__1_contigs__length_160647932,length=160647933>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_6__1_contigs__length_78844258,length=78844259>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_13__3_contigs__length_82050896,length=82050899>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_3__2_contigs__length_137956696,length=137956698>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_2__1_contigs__length_157594471,length=157594472>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_12__1_contigs__length_47609450,length=47609451>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_11__2_contigs__length_80009992,length=80009994>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_8__1_contigs__length_71271319,length=71271320>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_7__1_contigs__length_75018798,length=75018799>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_5__1_contigs__length_83128659,length=83128660>
##contig=<ID=Hap2_t_crist_hwy154_cen4119#2#Scaffold_10__2_contigs__length_75648701,length=75648703>
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=MAD,Number=1,Type=Integer,Description="Minimum site allele depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled probability estimate of the called genotype">
##FORMAT=<ID=GP,Number=1,Type=Float,Description="Genotype Probability, the log-scaled posterior probability of the called genotype">
##FORMAT=<ID=XD,Number=1,Type=Float,Description="eXpected Depth, background coverage as used for the Poisson model">
##FILTER=<ID=lowad,Description="Variant does not meet minimum allele read support threshold of 1">
##FILTER=<ID=lowdepth,Description="Variant has read depth less than 4">
##SAMPLE=<ID=2013FHA001>
##contig=<ID=Scaffold_10__2_contigs__length_75648701>
##contig=<ID=Scaffold_11__2_contigs__length_80009992>
##contig=<ID=Scaffold_12__1_contigs__length_47609450>
##contig=<ID=Scaffold_13__3_contigs__length_82050896>
##contig=<ID=Scaffold_1__1_contigs__length_160647932>
##contig=<ID=Scaffold_2__1_contigs__length_157594471>
##contig=<ID=Scaffold_3__2_contigs__length_137956696>
##contig=<ID=Scaffold_4__1_contigs__length_97222829>
##contig=<ID=Scaffold_5__1_contigs__length_83128659>
##contig=<ID=Scaffold_6__1_contigs__length_78844258>
##contig=<ID=Scaffold_7__1_contigs__length_75018798>
##contig=<ID=Scaffold_8__1_contigs__length_71271319>
##contig=<ID=Scaffold_9__2_contigs__length_79556474>
##bcftools_viewVersion=1.23-3-g34a49760+htslib-1.23-9-gacc28ac1
##bcftools_viewCommand=view -h 2013FHA001.vcf.gz; Date=Wed Mar  4 19:32:45 2026
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2013FHA001
```
bcftools +counts for same vcf:
```
bcftools +counts 2013FHA001.vcf.gz
[W::bcf_hdr_check_sanity] GP should be declared as Number=G
Number of samples: 1
Number of SNPs:    29859
Number of INDELs:  239728
Number of MNPs:    36580
Number of others:  78012
Number of sites:   413902
```

then can use bcftools merge to merge all vcfs? Do I need to rpoject this back onto snarls somehow to identify inversions? when I look at the snarls.json, there are 45265393 lines. Of those, there are 11779457 instances of the word "backward" which I think is related to the inversions because when I search one of the nodes of the inversion (9179874) it pops up in a set of parenthesis including backwards=true.
to merge vcfs:
```
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=mergevcf
#SBATCH --qos gompert-grn
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/merge-%A_%a.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/merge-%A_%a.out
#SBATCH --mem=200G

### load modules ###
module load bcftools
module load plink/2.0 

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/vg_vcf

### merge and index vcfs ###
bcftools merge *.vcf.gz -O z -threads 12 -o 2013FHA_merged.vcf.gz
bcftools index 2013FHA_merged.vcf.gz

### check missingness ###
plink --vcf 2013FHA_merged.vcf.gz --missing --out 2013FHA_merged_missingness_report
```
## GBS Data Alignment and Variant Calling from Pangenome with standard methods
First index the pangenome:

```
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

```
Then run mapping:
```
#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/bwamem_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/bwamem_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=bwamem
#SBATCH --qos gompert-grn
#SBATCH --array=0-1   # Job array when n is number of unique samples

### LOAD MODULES ###
#For this step, bwa and needed
module load bwa
module load samtools
module load bzip2

echo "Start Job"
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

### ASSIGN VARIABLES  ###
P=$(find /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/data -type l | sort | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
SAMPLE=$(basename $P | cut -d "." -f 1)
pangenome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.sv.gfa.fa.gz"

echo "P=$P"
echo "SAMPLE=$SAMPLE"
echo "pangenome=$pangenome"

### SET TMPDIR ###
WORKDIR="/scratch/general/nfs1/u6071015/timemaGBS/"
cd "$WORKDIR" 

### MAPPING ###
echo "Mapping ${SAMPLE}"

bwa mem -t 12 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" "$pangenome" <(bzip2 -dc "$P") > "${SAMPLE}.sam"
echo "Mapping complete for ${SAMPLE}"

### SAMTOOLS SORT AND FILTER ###
#convert sam to bam
samtools fixmate -O bam "${SAMPLE}.sam" "${SAMPLE}.bam"

samtools sort -@ 12 -o "${SAMPLE}.sorted.bam" "${SAMPLE}.bam"
samtools index "${SAMPLE}.sorted.bam"
samtools view -b -q 20 "${SAMPLE}.sorted.bam" > "${SAMPLE}.sorted.unique.bam"
samtools index "${SAMPLE}.sorted.unique.bam"

echo "samtools done for ${SAMPLE}"

```

and now variant calling:
```
#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/varcall_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/varcall_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=varcall
#SBATCH --qos gompert-grn

### LOAD MODULES ###
module load samtools
module load bcftools


### ASSIGN VARIABLES ###
BAMDIR="/scratch/general/nfs1/u6071015/timemaGBS/"
BAM_FILES=($(find "$BAMDIR" -type f -name "*.sorted.unique.bam" | sort ))
pangenome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.sv.gfa.fa.gz"
WORKDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf"
THREADS=20
MERGED="${WORKDIR}/FHA_all.unique.bam"
SORTED="${WORKDIR}/FHA_all.sorted.unique.bam"
OUTVCF="FHA_all.vcf"

### MERGE ALL BAMS FOR VARIANT CALLING ###
echo "Merge BAM files"
cd "$WORKDIR"

samtools merge -r -c -p -@ ${THREADS} "$MERGED" "${BAM_FILES[@]}"
samtools sort -@ 12 -o "$SORTED" "$MERGED"
samtools index "$SORTED"

echo "BAM files merged"

### VARIANT CALLING ###
echo "start variant calling"
cd "$WORKDIR"

bcftools mpileup -Ou --C 50 -d 500 -r "$CHR" -a DP,AD,ADF,ADR -q 20 -Q 30 -f "$pangenome" "$SORTED" | bcftools call -v -c -p 0.01 -P 0.001 -Ov -o "$OUTVCF"

echo "finished variant calling"
```

### Variant Filtering
Now for vcf filtering using Zachs filtering scripts ```vcfFilter.pl``` and ```filterSomeMoreL.pl``` from ChumashWGS(https://github.com/zgompert/ChumashWGSmapping), edited for these samples. 
Here is the script to run all:

```
#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/filter_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/filter_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=filter
#SBATCH --qos gompert-grn

### LOAD MODULES ###
module load R
module load bcftools

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf/

#plot stats on raw vcf
bcftools query \
       -f '%CHROM\t%POS\t%QUAL\t%INFO/MQBZ\t%INFO/SCBZ\t%INFO/RPBZ\t%INFO/DP\n' \
      FHA_all.vcf > FHA_all_site_metrics.txt
Rscript plot_histograms

#filtered the SNP set for coverage, missing data, and various tests of bias 
perl vcfFilter.pl FHA_all.vcf

#extracted the read depth per SNP and individual from the filtered vcf files
bcftools query -f '[%DP\t]\n' filtered0.5x_FHA_all.vcf | sed 's/\t$//' > depth.txt

#compute depth per individual and SNP to identify SNPs and individuals to drop
Rscript CovFilt.R

#Drop those individuals
perl filterSomeMoreL.pl filtered0.5x_FHA_all.vcf

#convert to genotype likelhood file
perl vcf2gl.pl 0.0 morefilter_0.5x_FHA_all.vcf

#obtain maximum likelihood of AF using expectation-maximization algorithm written by Zach, estpEM
estpEM -i FHA_all.gl -o FHA_all_0.5x_estpEM.txt -e 0.001 -m 50 -h 1

```
and the scripts within that script. here is vcfFilter.pl:
```
#!/usr/bin/perl

use warnings;
use strict;

# this program filters a vcf file based on overall sequence coverage, number of non-reference reads, number of alleles, and reverse orientation reads
# usage vcfFilter.pl infile.vcf
# change the marked variables below to adjust settings

#### stringency variables, edits as desired
## 602 inds, 2x
my $minCoverage = 1204; # minimum number of sequences; DP
my $minAltRds = 10; # minimum number of sequences with the alternative allele; AC
my $notFixed = 1.0; # removes loci fixed for alt; AF
my $bqrs = 3; # Z-score base quality rank sum test; BaseQRankSum
my $mqrs = 3; # Z-score mapping quality rank sum test; MQRankSum
my $rprs = 3; # Z-score read position rank sum test; ReadPosRankSum
my $mq = 30; # minimum mapping quality; MQ
my $miss = 60; # maximum number of individuals with no data = 10%
##### this set is for GBS
my $d;

my @line;

my $in = shift(@ARGV);
open (IN, $in) or die "Could not read the infile = $in\n";
$in =~ m/^([a-zA-Z_0-9\-]+)\.vcf$/ or die "Failed to match the variant file\n";
open (OUT, "> filtered2x_$1.vcf") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
	chomp;
	$flag = 1;
	if (m/^\#/){ ## header row, always write
		$flag = 1;
	}
	elsif (m/^\S+/){ ## this is a sequence line, you migh need to edit this reg. expr. used to be (m/^Sc/){
		$flag = 1;
		$d = () = (m/\d\/\d:0,0,0:0/g); ## for bcftools call
		if ($d >= $miss){
			$flag = 0;
			##print "fail missing : ";
		}
		if (m/[ACTGN]\,[ACTGN]/){ ## two alternative alleles identified
			$flag = 0;
			#print "fail allele : ";
		}
		@line = split(/\s+/,$_);
		if(length($line[3]) > 1 or length($line[4]) > 1){
			$flag = 0;
			#print "fail INDEL : ";
		}
		m/DP=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 < $minCoverage){
			$flag = 0;
			#print "fail DP : ";
		}
## bcftools call version
	
		m/DP4=\d+,\d+,(\d+),(\d+)/ or die "Syntax error DP4 not found\n";
		if(($1 + $2) < $minAltRds){
			$flag = 0;
		}
		m/AF1*=([0-9\.e\-]+)/ or die "Syntax error, AF not found\n";
		if ($1 == $notFixed){
			$flag = 0;
		#	print "fail AF : ";
		}

## bcftools call verions, these are p-values, use 0.01
		if(m/BQBZ=([0-9e\-\.]*)/){
			if (abs($1) > $bqrs){
				$flag = 0;
#				print "fail BQRS : ";
			}
		}
		if(m/MQBZ=([0-9e\-\.]*)/){
			if (abs($1) > $mqrs){
				$flag = 0;
#				print "fail MQRS : ";
			}
		}
		if(m/RPBZ=([0-9e\-\.]*)/){
			if (abs($1) > $rprs){
				$flag = 0;
#				print "fail RPRS : ";
			}
		}
		if(m/MQ=([0-9\.]+)/){
			if ($1 < $mq){
				$flag = 0;
#				print "fail MQ : ";
			}
		}
		else{
			$flag = 0;
			print "faile no MQ : ";
		}
		if ($flag == 1){
			$cnt++; ## this is a good SNV
		}
	}
	else{
		print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
		$flag = 0;
	}
	if ($flag == 1){
		print OUT "$_\n";
	}
}
close (IN);
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
```
Then I will extract read depth per individual and SNP, and idenitfy which ones to drop with CovFilt.R:

```
## compute depth per individual and SNP to identify SNPs and individuals to drop
## the idea is to get rid of low coverage individuals
## and SNPs with either very high coverage (3SD > mean) or very high variance in coverage
## across individuals

library(data.table)
d<-as.matrix(fread("depth.txt",header=FALSE))

## mean and SD by SNP
mnc<-apply(d,1,mean)
sdc<-apply(d,1,sd)

## CV 
cvc<-sdc/mnc
meancvc<-mean(cvc)
quantcvc<-quantile(cvc,probs=c(.5,.9,.95,.99,.999,1))

meanmnc3sd<-mean(mnc)+3*sd(mnc)
mean(mnc)
quantmnc<-quantile(mnc,probs=c(.5,.9,.95,.99,.999,1))

## for SNPs, keep if CV < 1.5 and mean < 21
keepSNPs<-as.numeric(mnc < 21 & cvc < 1.5)
keepSNPs_mean<-mean(keepSNPs)
keepSNPs_sum<-sum(keepSNPs)

## for individuals

mni<-apply(d,2,mean)
plot(sort(mni))
summary(mni)
quantile(mni,probs=c(.025,.99))

keepInds<-as.numeric(mni > 4 & mni < 40)
mean(keepInds)
sum(keepInds)

cat(sprintf("meancvc=%.4f\n", meancvc))
cat(sprintf("quantcvc=%s\n", paste(names(quantcvc), round(quantcvc, 4), sep="=", collapse=", ")))
cat(sprintf("meanmnc3sd=%.4f\n", meanmnc3sd))
cat(sprintf("meanmnc=%.4f\n", mean(mnc)))
cat(sprintf("quantmnc=%s\n", paste(names(quantile(mnc, probs=c(.5,.9,.95,.99,.999,1))), 
                                    round(quantile(mnc, probs=c(.5,.9,.95,.99,.999,1)), 4), 
                                    sep="=", collapse=", ")))
cat(sprintf("keepSNPs_mean=%.4f\n", mean(keepSNPs)))
cat(sprintf("keepSNPs_sum=%d\n", sum(keepSNPs)))
cat(sprintf("meanmni=%s\n", paste(names(summary(mni)), round(summary(mni), 4), sep="=", collapse=", ")))
cat(sprintf("quantmni=%s\n", paste(names(quantile(mni, probs=c(.025,.99))), 
                                    round(quantile(mni, probs=c(.025,.99)), 4), 
                                    sep="=", collapse=", ")))

write.table(file="KeepInds.txt",keepInds,row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(file="KeepSNPs.txt",keepSNPs,row.names=FALSE,col.names=FALSE,quote=FALSE)
```
and here is filterSomeMoreL.pl, based on the output of CovFilt.R, which resulted in vectors of 0s and 1s for SNPs and individuals to drop. For SNPs, I used flagged SNPs to remove with a CV > (around the 99.9th percentile) and mean coverage > 3 SDs above the mean. I flagged individuals with mean coverage < (2.5th percentile) or > (a bit above the 90th percentile). Drop the SNPs first; then convert to gl format. I used the following to drop the SNPs, resulting in the morefilter* vcf files.

```
#!/usr/bin/perl
# filter vcf files based on coverage

open(IN,"KeepSNPs.txt") or die "failed initial read\n";
while(<IN>){
	chomp;
	push(@keep,$_);
}
close(IN);


foreach $in (@ARGV){
	open (IN, $in) or die "Could not read the infile = $in\n";
	#$in =~ m/^([a-zA-Z0-9_]+\.vcf)$/ or die "Failed to match the variant file\n";
	open (OUT, "> morefilter_$1") or die "Could not write the outfile\n";


	while (<IN>){
		chomp;
		if (m/^\#/){ ## header row, always write
			$flag = 1;
		}
		elsif (m/^S+/){ ## this is a sequence line, you migh need to edit this reg. expr.
			$flag = shift(@keep);
			if ($flag == 1){
				$cnt++; ## this is a good SNV
			}
		}
		else{
			print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
			$flag = 0;
		}
		if ($flag == 1){
			print OUT "$_\n";
		}
	}
	close (IN);
	close (OUT);

	print "Finished filtering $in\nRetained $cnt variable loci\n";
}
```
Then run vcf2gl.pl. The output of this (FHA_all.gl) did not have the correct number of loci in the top line of the header (should be nind nloc), so I just manually edited the .gl file to include the correct number and then ran genotype likelihood estimation with Zach program estpEM(v1). Download main_estpEM.C, func_estpEM.C, and estpEM.H to folder and compile using this code:

```
g++ main_estpEM.C func_estpEM.C -lgsl -lgslcblas -lm -o estpEM
```
when running estpEM, there are a couple parameters that can be changed. 
-e 0.001 — convergence tolerance for the EM algorithm. It stops iterating when the allele frequency change between iterations is less than 0.001. This is actually already the default, so you'd only change it if you wanted stricter (e.g. 0.0001) or looser convergence.
-m 50 — maximum number of EM iterations before giving up. The default is 20, so this increases it to 50, giving the algorithm more chances to converge at tricky loci.
-h 2 — tells the program that your input file has 2 extra header lines to skip (beyond the required first line that specifies the number of individuals and loci). The default is 0.

Next steps from Zach's code:
Then obtain Bayesian estiamtes of genotypes using the allele frequency priors (under assuming HW genotype frequencies as prior expectations) for both the posterior mode and mean to compare. Compile the new C programs for the empirical Bayes genotype esimates, gl2genest.c and gl2genestMax.c. 
gcc gl2genest.c -lm -o gl2genest
gcc gl2genestMax.c -lm -o gl2genestMax


```
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
./gl2genestMax FHA_all_0.5x_estpEM.txt FHA_all.gl
## posterior mean
./gl2genest FHA_all_0.5x_estpEM.txt FHA_all.gl

```
The output files are cpntest_FHA_all.txt for the posterior mean and mlpntest_FHA_all.txt for the posterior mode. These contain 602 individuals (columns) and 5,340 SNPS (rows). 

### Local PCAs to identify inversions
Now for local PCAS using lostruct (https://github.com/petrelharp/local_pca?tab=readme-ov-file). Takes the output mean genotype likelihood file from last step

lostruct.R:
```
library(data.table)
devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)

# read the cpntest output (loci x individuals)
coded <- as.matrix(read.table("cpntest_FHA_all.txt"))
eigenstuff <- eigen_windows(coded, win=100, k=2)
windist <- pc_dist(eigenstuff, npc=2 )
fit2d <- cmdscale(windist, eig=TRUE, k=2 )
plot( fit2d$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(windist)) )

```
run_lostruct.sh
```
#!/bin/bash
#SBATCH --output=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/lostruct_%A_%a.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/logs/lostruct_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=lostruct
#SBATCH --qos gompert-grn

module load R

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf

Rscript lostruct.R

```
 
## Comparison across methods

To compare the success of calling across methods, we can use sveval (https://github.com/jmonlong/sveval) with vcfs from each method, or Zhang et al. 2025 then use survivor (https://www.github.com/fritzsedlazeck/SURVIVOR; version 1.0.3) (Jeffares et al., 2017) to identify homologous SV. Here is a survivor tutorial:
https://evomics.org/learning/population-and-speciation-genomics/2022-population-and-speciation-genomics/detecting-structural-variants-lab/

```
ls *vcf > sample_files
./SURVIVOR merge sample_files 1000 2 1 1 0 50 sample_merged.vcf
#maximum allowed distance of 1kb, supported by 2 callers, agree on the type (1) and on the strand (1) of the SV, at least 50bp
```
