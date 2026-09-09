# Timema-SV-Methods
This is code for a manuscript which investigates how different methods of structural variant (SV) detection and characterization vary, focused on inversions. Specifically, does SV size, frequency, or density have an effect on SV detection in different methods? Using whole-genome assemblies and genotyping-by-sequencing data from Gompert et al. 2025 (https://www.science.org/doi/10.1126/science.adp3745). Here is some helpful background on pangenomics: https://pangenome.github.io/

We compare three different approaches, all which have their own subfolder in this repository:
1) pangenome
2) comparative_alignment
3) local_pca

Overlap across the three methods is calculated and visualized in inversion.R in the main repository.

## PANGENOME

Starting off with making the pangenome using eight genomes (4 hwy154 and 4 refugio) in cactus minigraph (paper:https://www.nature.com/articles/s41587-023-01793-w, documentation: https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).

|  POP  | STRIPE |   ID  | SHORT |HAP|    SCIENCE BOUNDS    |
|-------|--------|-------|-------|---|----------------------|
|HWY154 |striped |cen4119| H GS  | 1 |24,457,103–39,030,359 |
|HWY154 |  green |cen4280| H GUS | 2 |24,803,527–44,121,870 |
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

Make the HWY154_REF.txt input file for Minigraph-cactus:

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

Cactus minigraph does require designating a reference, and you can actually designate multiple to use as the basis for VCF files (with command --vcfReference). The cactus publication runs on two different references and chooses the longest one. Other commands that I haven't run yet but could include if the graph seems weird, are --permissiveContigFilter and --noSplit (which disables chromosome splitting). Also, --vcfbub by default flattens the vcf and removes nested variants. If I want that not to happen, I have to set --vcfbub 0. This might make downstream annotation harder. Other things to note on this graph: minigraph only uses SVs > 50 bp in graph construction, and also clips out stretches of sequences >= 10 kb that do not align to minigraph. Now starting the minigraph pipeline with script ```run_cactus-pangenome.sh```:

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

Got message from job at time 11-10-2025 17:18:07: Job used more disk than requested. For CWL, consider increasing the outdirMin requirement, otherwise, consider increasing the disk requirement. Job 'unzip_gz' kind-unzip_gz/instance-b62o1k6q v1 used 102.02% disk (1.9 GiB [2029998080B] used, 1.9 GiB [1989727325B] requested).

Finished after about 10 hours in interactive job wiht 200G.


### Investigating Cactus Pangenome Output

halStats output for pangenomes with Refugio included:
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

From Science paper (Gompert et al. 2025), length of chromosomes:
*NOTE: The Refugio scaffolds were incorrect in this paper, I corrected them using pairwise synteny and correct ones are listed below*

Table S1: Homologous relationships among chromosome-size scaffolds for the T. cristinae
genomes. Chromosome 13 is the X sex chromosome. Abbreviations are as follows: Chr = chromosome,
U GS = unphased stripe genome from (82), R = Refugio, H = Hwy154, GS1 = striped
haplotype 1, GS2 = striped haplotype 2, GUS1 = green haplotype 1, and GUS2 = green haplotype.
4119, 4122=striped and 4120,4280=green. *reference for these analyses is H GS2*

|Chr |U GS |R GS1 |R GS2 |R GUS1 |R GUS2 |H GS1 |H GS2 |H GUS1 |H GUS2|
|----|-----|------|------|-------|------ |------|------|-------|------|
|  1 | 8483|   6  |   8  |  11   |   9   |  13  |  13  |  22   |  15  |
| 2  |14640|   8  |   5  |   7   |   5   |   5  |   6  |  23   |   1  |
| 3  |42935|   2  |   1  |   1   |   1   |   3  |   2  |  16   |   3  |
| 4  |42912|   1  |   1  |   1   |   1   |   1  |   1  |  64   |  35  |
| 5  |18722|  13  |  12  |  12   |  12   |  12  |  12  |   5   |  10  |
| 6  |9928 |   7  |   2  |   6   |   4   |   4  |   5  |  11   |  44  |
| 7  |10660|  10  |   7  |  10   |  11   |  10  |   8  |  54   |   7  |
| 8  |7748 |  11  |   9  |   3   |   3   |  11  |   4  |   7   |  23  |
| 9  |16151|   4  |   3  |   8   |  10   |   8  |   9  |  46   |  21  |
| 10 |14160|  12  |  10  |   4   |   6   |   7  |   7  |  15   |  16  |
| 11 |12033|   9  |   6  |   5   |   7   |   9  |  10  |   2   |  12  |
| 12 |12380|   5  |   4  |   9   |   8   |   6  |  11  |   1   |  36  |
| 13 |14101|   3  |  11  |   2   |   2   |   2  |   3  |  36   |   8  |

### Evaluating Pangenome quality with GRETL

Gretl github (https://github.com/MoinSebi/gretl) and manuscript (https://academic.oup.com/bioinformatics/article/41/1/btae755/7932228)

create conda environment:

```
cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome
conda create -n gretl_env
conda activate gretl_env
mamba install -c conda-forge -c bioconda gretl
#test
cargo test
```

run gretl:
```
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos=gompert-grn
#SBATCH --job-name=gretl
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/gretl-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/gretl-%j.out

module load miniforge3
conda activate gretl_env

SCAFF="Scaffold_4__1_contigs__length_97222829"
INPUT_DIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments"
WORK_DIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/gretl"

echo "Graph-based and hybrid stats"
./gretl stats -g ${INPUT_DIR}/${SCAFF}.gfa --pansn -o ${INPUT_DIR}/gretl_stats_${SCAFF}.txt

# Path-based statistics
#./gretl stats -g ${INPUT_DIR}/${SCAFF}.gfa --pansn -o ${INPUT_DIR}/gretl_pathstats_${SCAFF}.txt -p

echo "job done"
```

###  Pangenome Visualization with sequenceTubeMap 

Downloaded to local computer using these instructions: https://github.com/vgteam/sequenceTubeMap?tab=readme-ov-file
Then to run tube map on local computer in terminal:

```
cd ~/Desktop/GitHub/sequenceTubeMap
nvm use
npm run serve
```
I uploaded the .gbz and .gaf file from the cluster into the folder ~/Desktop/GitHub/sequenceTubeMap/exampleData following these instructions: https://github.com/vgteam/sequenceTubeMap/blob/master/doc/data.md. I went here to visualize: http://localhost:3000.

The .gbz gets mounted as the graph and haplotype, while the .gaf can be mounted as the reads. Make sure to index the .gaf file with tabix (htslib), and upload the index file to /exampleData/ as well. 

###  Pangenome Visualization with ODGI

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

### Pangenome Visualization with vg

Can run this in an interactive job to see specific nodes from inversions:

```
module load cactus/3.0.1

# try to make simpler graphs with just inversion nodes
vg find -x Scaffold_3__2_contigs__length_137956696.xg -n 6232384 -n 6232476 -c 3 | vg view -dp - | dot -Tsvg -o Scaffold3_subgraph6232384_6232476.svg
vg find -x Scaffold_13__3_contigs__length_82050896.xg -n 2993340 -n 2993702 -c 3 | vg view -dp - | dot -Tsvg -o Scaffold13_subgraph2993340_299370.svg
vg find -x Scaffold_13__3_contigs__length_82050896.xg -n 2993340 -n 2993702 | vg view -dp - | dot -Tsvg -o Scaffold13_subgraph2993340_299370_noc.svg
vg find -x Scaffold_4__1_contigs__length_97222829.xg -n 5937107 -n 6725313 -n 9179874 | vg view -dp - | dot -Tsvg -o Scaffold4_subgraph5937107_9179874_6725313.svg
vg find -x Scaffold_4__1_contigs__length_97222829.xg -n 5937107 -n 6725313 -n 9179874 | vg view -dp - | dot -Tsvg -o Scaffold4_subgraphall.svg
vg find -x Scaffold_4__1_contigs__length_97222829.xg -n 5937107 -n 5936984 -n 9179874 -n 11242967 -c 3| vg view -dn - -u | dot -Tsvg -o Scaffold4_subgraph5937107_9179874_5936984_11242967_c3.svg

vg find -x Scaffold_4__1_contigs__length_97222829.xg \
  $(awk '{printf "-n %s ", $1}' scaff4nodes.txt) -c 3 \
  | vg view -dn - -u | dot -Tsvg -o Scaffold4_subgraphall_c3.svg

vg find -x Scaffold_9__2_contigs__length_79556474.xg \
  $(awk '{printf "-n %s ", $1}' scaff9nodes.txt) -c 3 \
  | vg view -dn - -u | dot -Tsvg -o Scaffold9_subgraphall_c3.svg
```

### Using Pantree to describe SVs in the pangenome
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
Now to summarize the output vcf using code from their paper (https://github.com/ShenghanZhang1123/graph_var_analysis/blob/main/notebooks/generating_data_analysis.ipynb) in an adapted script I wrote called pantree_summary.py:

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

mkdir summary

#summarize SVs (edited from pantree manuscript, puts output in /summary subdir of working directory)
python pantree_summary.py --vcf /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/${SCAFF}_pantree.vcf.gz --chrom ${SCAFF}

```

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
### Projected pantree output back into genome coordinate spaces

Script to search all inversion vcfs for genome coordinates corresponding to inversion nodes. Some genomes do not have positions for the nodes, this means that their path does not pass through that node.
```
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=summarizecoords
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizecoords%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizecoords-%j.out

module load cactus/3.0.1

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree

# Output table header
echo -e "scaffold\tvariant_id\tgenome\tstart_node\tstart_pos\tend_node\tend_pos" > inversion_coordinates.tsv

# Loop over each scaffold's inversion VCF
for vcf in *_pantree_inversions_only.vcf.gz; do

    # Extract scaffold name/number from filename
    scaffold=$(basename "$vcf" .vcf.gz | sed 's/_pantree_inversions_only//')
    
    # Path to corresponding xg index, vg file, and display paths file
    xg="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/${scaffold}.xg"
    vg_file="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/${scaffold}.vg"
    paths_file="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/chrom-alignments/display_paths_${scaffold}.txt"
    
    # Check vg file exists
    if [[ ! -f "$vg_file" ]]; then
        echo "WARNING: $vg_file not found, skipping $scaffold" >&2
        continue
    fi

    # Build xg index if it doesn't exist
    if [[ ! -f "$xg" ]]; then
        echo "INFO: $xg not found, building from $vg_file..." >&2
        vg index -x "$xg" "$vg_file"
        if [[ $? -ne 0 ]]; then
            echo "ERROR: failed to build $xg, skipping $scaffold" >&2
            continue
        fi
        echo "INFO: $xg built successfully" >&2
    else
        echo "INFO: $xg already exists, skipping index build" >&2
    fi

    # Check display paths file exists
    if [[ ! -f "$paths_file" ]]; then
        echo "WARNING: $paths_file not found, skipping $scaffold" >&2
        continue
    fi

    # Loop over each inversion in the VCF (skip header lines, decompress on the fly)
    while IFS=$'\t' read -r chrom pos id ref alt qual filter info format samples; do
        
        # Extract start and end nodes from the ID field e.g. >2993340<2993702
        start_node=$(echo "$id" | grep -oP '(?<=[><])\d+' | head -1)
        end_node=$(echo "$id"   | grep -oP '(?<=[><])\d+' | tail -1)
        
        # Loop over each genome path
        while IFS= read -r genome_path; do
            
            # Skip empty lines
            [[ -z "$genome_path" ]] && continue
            
            # Query start node
            start_result=$(vg find -x "$xg" -n "$start_node" -P "$genome_path" 2>/dev/null)
            start_pos=$(echo "$start_result" | awk '{print $2}')
            
            # Query end node
            end_result=$(vg find -x "$xg" -n "$end_node" -P "$genome_path" 2>/dev/null)
            end_pos=$(echo "$end_result" | awk '{print $2}')
            
            # If either position is empty, mark as missing
            [[ -z "$start_pos" ]] && start_pos="."
            [[ -z "$end_pos" ]]   && end_pos="."
            
            # Write to output table
            echo -e "${scaffold}\t${id}\t${genome_path}\t${start_node}\t${start_pos}\t${end_node}\t${end_pos}" \
                >> inversion_coordinates.tsv

        done < "$paths_file"

    done < <(zcat "$vcf" | grep -v "^#")

done

echo "Done. Output in inversion_coordinates.tsv"

```

And then to get the allele length and genotype information from each of the vcfs:

```
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos gompert-grn
#SBATCH --job-name=summarizelengths
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizelengthsgenos%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/logs/summarizelengthsgenos-%j.out

module load bcftools

FIRST=1

for vcf in *_pantree_inversions_only.vcf.gz; do
    scaffold=$(basename "$vcf" .vcf.gz | sed 's/_pantree_inversions_only//')

    SAMPLES=$(bcftools query -l "$vcf" \
      | grep 't_crist' \
      | grep -v 'MINIGRAPH' \
      | tr '\n' ',' \
      | sed 's/,$//')

    if [[ -z "$SAMPLES" ]]; then
        echo "WARNING: no t_crist samples found in $vcf, skipping" >&2
        continue
    fi

    if [[ $FIRST -eq 1 ]]; then
        # Print header on first VCF only
        bcftools query -s "$SAMPLES" -H \
          -f '%ID\t%INFO/RC\t%INFO/AC\t%INFO/TP\t%INFO/NIA\t%INFO/AN[\t%GT]\n' \
          "$vcf"
        FIRST=0
    else
        # No header for subsequent VCFs
        bcftools query -s "$SAMPLES" \
          -f '%ID\t%INFO/RC\t%INFO/AC\t%INFO/TP\t%INFO/NIA\t%INFO/AN[\t%GT]\n' \
          "$vcf"
    fi
done > all_scaffolds_inversions_tcrist_genotypes_allNR.tsv

awk 'BEGIN{OFS="\t"}
  NR==1 {
    # Fix header: remove NR col (8), rename NR_count col (9)
    for (i=1; i<=NF; i++) {
      if (i==8) continue
      if (i==9) printf "NR_bp"
      else printf "%s", $i
      if (i!=NF) printf "\t"
    }
    printf "\n"
    next
  }
  {
    nr_val = $8
    nr_bp = (nr_val == ".") ? 0 : length(nr_val)
    for (i=1; i<=NF; i++) {
      if (i==8) continue
      if (i==9) printf "%s", nr_bp
      else printf "%s", $i
      if (i!=NF) printf "\t"
    }
    printf "\n"
  }
' all_scaffolds_inversions_tcrist_genotypes_allNR.tsv > all_scaffolds_inversions_tcrist_genotypes_NRbp.tsv
```
### test inversion calling with INVPG_annot
program github:https://github.com/SandraLouise/INVPG_annot

to run:
```
invpg -v HWY154_REF_4119Hap2.vcf -g HWY154_REF_4119Hap2.gfa -o invpg_HWY154_REF_4119Hap2.vcf -k -m 0.5 -d 10
```
## LOCAL PCA

### GBS Data Alignment and Variant Calling with standard methods

I will be using both GBS data from HWY 154 and Refugio to quantify inversions using a local pca approach. This is the same data from the 2025 Science paper.

602 FHA individuals located:
/uufs/chpc.utah.edu/common/home/gompert-group3/data/sheffield/timema/2013fha_gwas/02_ids_reads/cristinae/2013*bz2

238 Refugio individuals located:
/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/clines/2016_gwas_trad_Patrik_clines/parsed/16*fastq

ln -s /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/clines/2016_gwas_trad_Patrik_clines/parsed/16*fastq /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/REF_data/

This was run with both the pangenome and the GHS2 reference. Joint variant calling did not work for the pangenome, so we ultimately just used the one reference.

Run mapping:
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
#SBATCH --array=1-602   # Job array when n is number of unique samples

### LOAD MODULES ###
#For this step, bwa and needed
module load bwa
module load samtools

echo "Start Job"
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

### ASSIGN VARIABLES  ###
P=$(find /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/data -type l | sort | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
SAMPLE=$(basename $P | cut -d "." -f 1)
pangenome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.sv.gfa.fa.gz"
genome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked"
echo "P=$P"
echo "SAMPLE=$SAMPLE"
echo "pangenome=$pangenome"
echo "genome=$genome"

### SET TMPDIR ###
WORKDIR="/scratch/general/nfs1/u6071015/timemaGBS/"
cd "$WORKDIR"

### MAPPING ###
echo "Mapping ${SAMPLE}"

#replicated Science paper
bwa aln -n 4 -k 2 -l 20 -q 10 "$genome" <(bzcat "$P") > "${SAMPLE}_aligned.sai"
bwa samse "$genome" "${SAMPLE}_aligned.sai"  <(bzcat "$P") | \
    samtools view -bS -q 1 - | \
    samtools sort - > "${SAMPLE}.sorted.unique.bam"
samtools index "${SAMPLE}.sorted.unique.bam"

echo "Mapping complete for ${SAMPLE}"


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
#SBATCH --mem=400G

### LOAD MODULES ###
module load samtools/1.16
module load bcftools/1.16
#ran previously with bcftools/1.23

### ASSIGN VARIABLES ###
BAMDIR="/scratch/general/nfs1/u6071015/timemaGBS/"
BAM_FILES=($(find "$BAMDIR" -type f -name "*.sorted.unique.bam" | sort ))
pangenome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus_pangenome/HWY154_REF_4119Hap2/HWY154_REF_4119Hap2.sv.gfa.fa.gz"
genome="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked"
WORKDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf"
THREADS=12
MERGED="${WORKDIR}/FHA_all_oneref.unique.bam"
SORTED="${WORKDIR}/FHA_all_oneref.sorted.unique.bam"
OUTVCF="FHA_all_oneref.vcf"

### MERGE ALL BAMS FOR VARIANT CALLING ###
echo "Merge BAM files"
cd "$WORKDIR"

#samtools merge -f -r -c -p -@ ${THREADS} "$MERGED" "${BAM_FILES[@]}"
#samtools sort -@ 12 -o "$SORTED" "$MERGED"
#samtools index "$SORTED"
#samtools flagstat "$SORTED"

echo "BAM files merged"

### VARIANT CALLING ###
echo "start variant calling"
cd "$WORKDIR"

#same commands as science paper
bcftools mpileup -Ou -d 500 -a DP,AD,ADF,ADR -Q 30 -q 20 --skip-indels -f "$genome" "${BAM_FILES[@]}" | \
    bcftools call -v -c -p 0.01 -P 0.001 -Ov -o "$OUTVCF"

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
Then I will extract read depth per individual and SNP, and idenitfy which ones to drop with CovFilt.R. First I need to edit the hardcoded filters in CovFilt.R based off of the metrics in the data:
```
Retained 335645 variable loci
[1] 6.106647
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
 0.006599  5.095792  6.151516  6.106647  7.078124 11.268322 
    2.5%      99% 
2.717672 9.889495 
[1] 0.9152824
[1] 551
meancvc=0.5862
quantcvc=50%=0.5623, 90%=0.7579, 95%=0.8055, 99%=0.8757, 99.9%=0.9918, 100%=1.5029
meanmnc3sd=14.6447
meanmnc=6.1066
quantmnc=50%=5.7292, 90%=9.4419, 95%=10.6578, 99%=13.5714, 99.9%=23.6611, 100%=419.5266
keepSNPs_mean=0.9985
keepSNPs_sum=335142
meanmni=Min.=0.0066, 1st Qu.=5.0958, Median=6.1515, Mean=6.1066, 3rd Qu.=7.0781, Max.=11.2683
quantmni=2.5%=2.7177, 99%=9.8895
Finished filtering filtered2x_FHA_all_oneref.vcf
Retained 335142 variable loci
```
CovFilt.R

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

## for SNPs, keep if CVC < 99.9% and mean < 3sdmeanmnc
keepSNPs<-as.numeric(mnc < 15 & cvc < 1)
keepSNPs_mean<-mean(keepSNPs)
keepSNPs_sum<-sum(keepSNPs)

## for individuals

mni<-apply(d,2,mean)
plot(sort(mni))
summary(mni)
quantile(mni,probs=c(.025,.99))

#keep if mena coverage (mni) is >2.5% and <90%
keepInds<-as.numeric(mni > 2.7 & mni < 9)
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
I ended up with full file paths in the sample names. need to edit:
```
module load bcftools
bcftools query -l  morefilter_2x_FHA_all_oneref_v2.vcf | xargs -I{} basename {} .sorted.unique.bam > new_samples.txt
bcftools reheader -s new_samples.txt -o morefilt_rehead_2x_FHA_all_oneref_v2.vcf morefilter_2x_FHA_all_oneref_v2.vcf
bcftools query -l morefilt_rehead_2x_FHA_all_oneref_v2.vcf
bcftools query -f '%CHROM\t%POS\n' morefilt_rehead_2x_FHA_all_oneref_v2.vcf | sed 's/.*|s//' > positions.txt

```
Then run vcf2gl.pl. The output of this (FHA_all.gl) did not have the correct number of loci in the top line of the header (should be nind nloc), so I just manually edited the .gl file to include the correct number (602 FHA_all.gl  )and then ran genotype likelihood estimation with Zach program estpEM(v1). Download main_estpEM.C, func_estpEM.C, and estpEM.H to folder and compile using this code:

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
./gl2genestMax FHA_all_v2_estpEM.txt  FHA_all.gl
## posterior mean
./gl2genest FHA_all_v2_estpEM.txt  FHA_all.gl

```
The output files are cpntest_FHA_all.txt for the posterior mean and mlpntest_FHA_all.txt for the posterior mode. These contain 602 individuals (columns) and 5,340 SNPS (rows). 

### Local PCAs to identify inversions with Lostruct
Now for local PCAS using lostruct (https://github.com/petrelharp/local_pca?tab=readme-ov-file). Takes the output mean genotype likelihood file from last step.
This also requires a positions.txt file, created in the sbatch script to run lostruct:

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

module load bcftools
module load R

SCRIPTDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf"
WORKDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/GBS/bcftools_vcf/FHA_alignedtocen4119hap2"

cd ${WORKDIR}

#make positions file
bcftools query -f '%CHROM\t%POS\n' morefilter_2x_REF_all_oneref.vcf | sed 's/.*|s//' > positions_raw.txt
awk -F'\t' 'BEGIN {
    map[1]=4; map[2]=3; map[3]=13; map[4]=8; map[5]=6;
    map[6]=2; map[7]=10; map[8]=7; map[9]=9; map[10]=11;
    map[11]=12; map[13]=1
}
{
    match($1, /Scaffold_([0-9]+)__/, arr)
    n = arr[1]+0
    $1 = (n in map) ? "Chr" map[n] : "UNMAPPED_" $1
    print $1 "\t" $2
}' positions_raw.txt > positions.txt


#Usage: Rscript ${SCRIPTDIR}/localpca_manyMDSaxes_v2.R <input_file> <output_prefix> <window_size_snps> <n_axes>
Rscript ${SCRIPTDIR}/localpca_manyMDSaxes_withfiltering.R cpntest_FHA_all.txt FHA_all 100 10

```
## COMPARATIVE ALIGNMENT

### Pairwise comparison in Progressive Cactus

We are also going to call SVs from the pairwise comparisons, specifically focusing on comparisons between the reference haplotype used for the pangenome (HWY154 Stripe Haplotype2) and the other haplotypes. For this, I am creating softlinks in ''/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus'' to existing hal files, and need to make the hal file for H154 Stripe 2/Refugio Stripe 1 pair. To do so, I run the script run_cactus.sh in the directory. The script also requires an input file with the genome name and path. For example, cactusTcrGSH2_TcrGSR1.txt:

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

### Whole genome assembly Pairwise comparison with syRI

#### install
```
conda activate syRI
conda install python=3.11
conda install cython numpy scipy pandas biopython psutil matplotlib plotsr
conda install -c conda-forge python-igraph 
conda install -c bioconda pysam 
conda install -c bioconda longestrunsubsequence
conda install -c bioconda syri

```

#### prepare fasta files
from documentation (https://schneebergerlab.github.io/syri/pipeline.html): Ideally, syri expects that the homologous chromosomes in the two genomes would have exactly same chromosome id. Therefore, it is recommended that the user pre-processes the fasta files to ensure that homologous chromosomes have exactly the same id in both fasta files corresponding to the two genomes. In case, that is not the case, syri would try to find homologous genomes using whole genome alignments, but that method is heuristical and can result in suboptimal results. Also, it is recommended that the two genomes (fasta files) should have same number of chromosomes.

python script to rename based on table from science:
```
#!/usr/bin/env python3
"""
Rename fasta scaffold headers to homologous chromosome IDs for SyRI.
Scaffold_N__... -> >Chr1, >Chr2, etc. based on lookup table.
"""

import os
import re

# --- Configuration ---
INPUT_DIR = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes"
OUTPUT_DIR = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed"

# Chromosome mapping table: Chr -> {genome: scaffold_number}
# Columns: Chr | UGS | RGS1 | RGS2 | RGUS1 | RGUS2 | HGS1 | HGS2 | HGUS1 | HGUS2
CHROM_TABLE = [
    #Chr   UGS    RGS1  RGS2  RGUS1  RGUS2  HGS1  HGS2  HGUS1  HGUS2
    ( 1,  8483,    12,   10,    4,     6,    13,   13,    22,    15),
    ( 2, 14640,     6,    8,   11,     9,     5,    6,    23,     1),
    ( 3, 42935,     2,    1,    1,     1,     3,    2,    16,     3),
    ( 4, 42912,     1,    1,    1,     1,     1,    1,    64,    35),
    ( 5, 18722,     7,    2,    6,     4,    12,   12,     5,    10),
    ( 6,  9928,     8,    5,    7,     5,     4,    5,    11,    44),
    ( 7, 10660,    10,    7,   10,    11,    10,    8,    54,     7),
    ( 8,  7748,    11,    9,    3,     3,    11,    4,     7,    23),
    ( 9, 16151,     5,    4,    9,     8,     8,    9,    46,    21),
    (10, 14160,     4,    3,    8,    10,     7,    7,    15,    16),
    (11, 12033,     9,    6,    5,     7,     9,   10,     2,    12),
    (12, 12380,    13,   12,   12,    12,     6,   11,     1,    36),
    (13, 14101,     3,   11,    2,     2,     2,    3,    36,     8),
]

# Genome file mapping: filename -> (genome_name, column index in CHROM_TABLE)
# hapltypes with scaffold 1 fusion excluded from these
GENOMES = {
    "t_crist_hwy154_cen4119_hap1.fasta.masked": ("HGS1",  6),
    "t_crist_hwy154_cen4119_hap2.fasta.masked": ("HGS2",  7),
    "t_crist_hwy154_cen4280_hap1.fasta.masked": ("HGUS1", 8),
    "t_crist_hwy154_cen4280_hap2.fasta.masked": ("HGUS2", 9),
    "t_crist_refug_cen4122_hap1.fasta.masked":  ("RGS1", 2),

}

os.makedirs(OUTPUT_DIR, exist_ok=True)

for filename, (genome_name, col_idx) in GENOMES.items():
    input_path = os.path.join(INPUT_DIR, filename)
    output_path = os.path.join(OUTPUT_DIR, filename)

    scaffold_to_chr = {}
    for row in CHROM_TABLE:
        chr_num = row[0]
        scaffold_num = row[col_idx]
        scaffold_to_chr[scaffold_num] = f"Chr{chr_num}"

    print(f"Processing {genome_name}: {filename}")
    print(f"  Scaffold -> Chr mapping: {scaffold_to_chr}")

    renamed = 0
    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                # Extract N from ">Scaffold_N__..."
                match = re.match(r">Scaffold_(\d+)__", line)
                if match:
                    scaffold_num = int(match.group(1))
                    if scaffold_num in scaffold_to_chr:
                        new_header = f">{scaffold_to_chr[scaffold_num]}\n"
                        fout.write(new_header)
                        renamed += 1
                    else:
                        print(f"  WARNING: Scaffold_{scaffold_num} not found in mapping table, keeping original header")
                        fout.write(line)
                else:
                    print(f"  WARNING: Could not parse header: {line.strip()}, keeping original")
                    fout.write(line)
            else:
                fout.write(line)

    print(f"  Done: {renamed}/13 chromosomes renamed -> {output_path}\n")

print("All genomes processed.")
print(f"Renamed fastas are in: {OUTPUT_DIR}")
```
To run it:
```
python3 rename_fasta_headers.py
```
then I want to fuse the chromosomes 3&4 in the reference to match the three fused haplotypes. I believe the orienation is like this, 3 different version os fhte reference GSH2:
RGUS1 version: Chr3 + Chr 4 fused
RGUS2 version: Chr4 + Chr3 inverted
RGS2: Chr4+Chr3 inverted

length chr 2 (3 in final; HGS2): 157594472 
length chr 1 (4 in final; (HGS2): 160647933 
length chr1 (RGS2): 318039304
length chr1 (RGUS1): 320397043
length chr1(RGUS2): 322377046

So then the code to flip these would be this python script ```python3 create_fused_references.py```:

```
#!/usr/bin/env python3
"""
Create three versions of the HGS2 reference genome with Chr3 and Chr4 fused,
matching the fusion arrangements in RGS2, RGUS1, and RGUS2.

Fusion arrangements:
  RGUS1 version: Chr3 + Chr4  (Chr3 first, Chr4 appended)
  RGUS2 version: Chr4 + Chr3  (Chr4 first, Chr3 reverse complemented appended)
  RGS2  version: Chr4 + revcomp(Chr3)  (Chr4 first, Chr3 reverse complemented appended)

The fused scaffold is named Chr3 in all versions (as Chr4 is absorbed into it),
and Chr4 is removed as a separate entry.
"""

import os

INPUT_FASTA = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_hwy154_cen4119_hap2.fasta.masked"
OUTPUT_DIR  = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed"

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]

def parse_fasta(filepath):
    """Read fasta into ordered list of (header, sequence) tuples."""
    records = []
    header = None
    seq_parts = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            records.append((header, "".join(seq_parts)))
    return records

def write_fasta(records, filepath, line_width=50):
    with open(filepath, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + "\n")
records = parse_fasta(INPUT_FASTA)

# Extract Chr3 and Chr4, keep everything else in order
chr3_seq = None
chr4_seq = None
other_records = []

for header, seq in records:
    name = header.split()[0]
    if name == "Chr3":
        chr3_seq = seq
        print(f"  Found Chr3: {len(seq):,} bp")
    elif name == "Chr4":
        chr4_seq = seq
        print(f"  Found Chr4: {len(seq):,} bp")
    else:
        other_records.append((header, seq))

if chr3_seq is None or chr4_seq is None:
    raise ValueError("Could not find Chr3 and/or Chr4 in input fasta. Check header names.")

# Define the three fusion versions
versions = {
    "RGUS1": (chr3_seq + chr4_seq,         "Chr3 + Chr4"),
    "RGUS2": (chr4_seq + chr3_seq,         "Chr4 + Chr3"),
    "RGS2":  (chr4_seq + revcomp(chr3_seq),"Chr4 + revcomp(Chr3)"),
}

for genome_name, (fused_seq, description) in versions.items():
    # Build record list: replace Chr3 with fused sequence, drop Chr4
    output_records = []
    for header, seq in records:
        name = header.split()[0]
        if name == "Chr3":
            output_records.append((f"Chr3", fused_seq))
            print(f"\n  {genome_name}: fusing {description} -> {len(fused_seq):,} bp")
        elif name == "Chr4":
            pass  # absorbed into Chr3
        else:
            output_records.append((header, seq))

    outfile = os.path.join(OUTPUT_DIR, f"t_crist_hwy154_cen4119_hap2.fasta.masked.fused_{genome_name}")
    write_fasta(output_records, outfile)

print("\nDone. Created 3 fused reference versions:")
```
make genomes file for plotsr:

```
# Create genomes file
echo "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_hwy154_cen4119_hap2.fasta.masked	HGS2
/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_hwy154_cen4119_hap1.fasta.masked.reoriented	HGS1" > genomes.txt
```

### Run SyRI pipeline
```
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --mem=100G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos=gompert-grn
#SBATCH --job-name=syri
#SBATCH -e /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/syri-%j.err
#SBATCH -o /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/logs/syri-%j.out

module load miniforge3
conda activate syRI

#set paths
cwd="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/syri" 
REFGENOME="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_hwy154_cen4119_hap2.fasta.masked"
QRYGENOME="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed/t_crist_hwy154_cen4119_hap1.fasta.masked.reoriented"
OUT="syri_TcrGSH2_TcrGSH1"

#perform whole genome alignment
minimap2 -ax asm5 --eqx ${REFGENOME} ${QRYGENOME} > ${OUT}.sam

#runSyRI, -k means keep intermediate files, -F is input format, .sam
syri -c ${OUT}.sam -r ${REFGENOME} -q ${QRYGENOME} -k -F S --nosnp 

#Plotting genomic structures predicted by SyRI
plotsr ${OUT}_syri.out ${REFGENOME} ${QRYGENOME} -H 8 -W 5
```

I got an error that large proportions of some of the genomes were inverted, so I wrote a python script to invert the chromosomes that were flagged. in this script, I have to change the name of the genome and the chromosomes I want inverted. it is located in /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/renamed, and I can run it using ```python3 invert_chroms.py ```.

### merge syri.out to evaluate shared inversions

I have output from syri pairwise alignments and sv calling (syri.out described here https://schneebergerlab.github.io/syri/fileformat.html) from seven different genomes aligned to the reference. They are located in their own folders: syri/syri_TsrGSH2_${QRY}/syri.out. I will write an R script to run with sbatch that loops through these files and finds all inversions (TYPE=INV, INVTR, INVDP):

R script ```group_syri_inversions.R```:

```
#!/usr/bin/env Rscript
# Group inversions from syri pairwise alignments across
# multiple query genomes vs. a common reference.
#
# Usage: Rscript analyze_syri_inversions.R
# Output directory: ./inversion_analysis/

## Load packages
library(dplyr)
library(purrr)

#### Config ####
# Query genome identifiers (adjust as needed)
QUERY_GENOMES <- c(
  "TcrHGS1",
  "TcrRGUS1",
  "TcrRGUS2",
  "TcrRGS2",
  "TcrRGS1",
  "TcrHGUS2",
  "TcrHGUS1"
)

# Base path pattern for syri output files
SYRI_PATH_TEMPLATE <- "syri_TcrGSH2_%s/syri.out"

# Inversion annotation types to include
INV_TYPES <- c("INV", "INVTR",  "INVDP")

# Output directory
OUT_DIR <- "inversion_analysis"

# Column names for syri.out 
SYRI_COLS <- c(
  "ref_chr", "ref_start", "ref_end",
  "ref_seq", "qry_seq",
  "qry_chr", "qry_start", "qry_end",
  "unique_id", "parent_id",
  "ann_type", "copy_status"
)

#### READ & FILTER INVERSIONS ####
read_inversions <- function(genome_id) {
  path <- sprintf(SYRI_PATH_TEMPLATE, genome_id)
  if (!file.exists(path)) {
    warning(sprintf("File not found for genome %s: %s", genome_id, path))
    return(NULL)
  }
  df <- tryCatch(
    read.table(path, sep = "\t", header = FALSE,
               col.names = SYRI_COLS, quote = "",
               comment.char = "#", fill = TRUE,
               stringsAsFactors = FALSE),
    error = function(e) {
      warning(sprintf("Error reading %s: %s", path, e$message))
      return(NULL)
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Keep only inversion-related annotations
  df <- df %>%
    filter(ann_type %in% INV_TYPES) %>%
    mutate(
      genome   = genome_id,
      ref_start = as.integer(ref_start),
      ref_end   = as.integer(ref_end),
      qry_start = as.integer(qry_start),
      qry_end   = as.integer(qry_end),
    ) %>%
    mutate(
      ref_start = pmin(ref_start, ref_end),
      ref_end   = pmax(ref_start, ref_end),
      inv_size  = ref_end - ref_start + 1L
    ) %>%
    select(genome, ref_chr, ref_start, ref_end,
           qry_chr, qry_start, qry_end,
           unique_id, parent_id, ann_type, inv_size)
  df
}

inv_list <- map(QUERY_GENOMES, read_inversions)
names(inv_list) <- QUERY_GENOMES
inv_list  <- compact(inv_list)   # drop NULLs
all_inv <- bind_rows(inv_list)
n_genomes <- length(inv_list)
message(sprintf("  Loaded %d inversion records across %d genomes.",
                nrow(all_inv), length(inv_list)))

size_sum <- all_inv %>%
  summarise(min    = min(inv_size),
            median = median(inv_size),
            mean   = mean(inv_size),
            max    = max(inv_size))
message(sprintf("  Size summary (bp): min=%d  median=%.0f  mean=%.0f  max=%d",
                size_sum$min, size_sum$median, size_sum$mean, size_sum$max))

# Save raw table
write.table(all_inv,
            file.path(OUT_DIR, "all_inversions_raw.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
```

all_inversions_raw.tsv is then used as the input file for ```inversion.R``` script where the inversions are grouped by reciprocal overlap and then analyzed and plotted further.

```
#!/usr/bin/env Rscript
# Group inversions from syri pairwise alignments across
# multiple query genomes vs. a common reference.
#
# Usage: Rscript analyze_syri_inversions.R
# Output directory: ./inversion_analysis/

## Load packages
library(dplyr)
library(purrr)

#### Config ####
# Query genome identifiers (adjust as needed)
QUERY_GENOMES <- c(
  "TcrHGS1",
  "TcrRGUS1",
  "TcrRGUS2",
  "TcrRGS2",
  "TcrRGS1",
  "TcrHGUS2",
  "TcrHGUS1"
)

# Base path pattern for syri output files
SYRI_PATH_TEMPLATE <- "syri_TcrGSH2_%s/syri.out"

# Inversion annotation types to include
SV_TYPES <- c("DUP")

# Output directory
OUT_DIR <- "inversion_analysis"

# Column names for syri.out 
SYRI_COLS <- c(
  "ref_chr", "ref_start", "ref_end",
  "ref_seq", "qry_seq",
  "qry_chr", "qry_start", "qry_end",
  "unique_id", "parent_id",
  "ann_type", "copy_status"
)

#### READ & FILTER INVERSIONS ####
read_SV <- function(genome_id) {
  path <- sprintf(SYRI_PATH_TEMPLATE, genome_id)
  if (!file.exists(path)) {
    warning(sprintf("File not found for genome %s: %s", genome_id, path))
    return(NULL)
  }
  df <- tryCatch(
    read.table(path, sep = "\t", header = FALSE,
               col.names = SYRI_COLS, quote = "",
               comment.char = "#", fill = TRUE,
               stringsAsFactors = FALSE),
    error = function(e) {
      warning(sprintf("Error reading %s: %s", path, e$message))
      return(NULL)
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Keep only SV-related annotations
  df <- df %>%
    filter(ann_type %in% SV_TYPES) %>%
    mutate(
      genome   = genome_id,
      ref_start = as.integer(ref_start),
      ref_end   = as.integer(ref_end),
      qry_start = as.integer(qry_start),
      qry_end   = as.integer(qry_end),
    ) %>%
    mutate(
      ref_start = pmin(ref_start, ref_end),
      ref_end   = pmax(ref_start, ref_end),
      inv_size  = ref_end - ref_start + 1L
    ) %>%
    select(genome, ref_chr, ref_start, ref_end,
           qry_chr, qry_start, qry_end,
           unique_id, parent_id, ann_type, inv_size)
  df
}

sv_list <- map(QUERY_GENOMES, read_SV)
names(sv_list) <- QUERY_GENOMES
sv_list  <- compact(sv_list)   # drop NULLs
all_sv <- bind_rows(sv_list)
n_genomes <- length(sv_list)
message(sprintf("  Loaded %d SV records across %d genomes.",
                nrow(all_sv), length(sv_list)))

size_sum <- all_sv %>%
  summarise(min    = min(inv_size),
            median = median(inv_size),
            mean   = mean(inv_size),
            max    = max(inv_size))
message(sprintf("  Size summary (bp): min=%d  median=%.0f  mean=%.0f  max=%d",
                size_sum$min, size_sum$median, size_sum$mean, size_sum$max))

# Save raw table
write.table(all_sv, file.path(OUT_DIR, "all_duplications_raw.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
```

## OVERLAP: Comparison across methods
Used custom r script ```inversions.R```
