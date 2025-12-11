# Timema-SV-Methods
Identifying how different methods of structural variant (SV) detection identify SVs, and whether Does SV size, frequency, or age have an effect on SV detection in different methods? Using data from this Gompert et al. 2025 (https://www.science.org/doi/10.1126/science.adp3745), and specifically focused on inversions and translocations. Here is some helpful background on pangenomics: https://pangenome.github.io/

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
Possible analysis can be done in odgi: https://odgi.readthedocs.io/en/latest/

### Investigating Cactus Output

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
4119, 4122=striped and 4120,4280=green

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

Downloaded using conda on kingspeak, because granite seems to not be working

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

WORKDIR="/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus/chrom-alignments"
SCAFFOLD="Scaffold_4__1_contigs__length_97222829"

cd $WORKDIR

#odgi build -g ${SCAFFOLD}.gfa  -o ${SCAFFOLD}.og -O

odgi sort -i ${SCAFFOLD}.og --threads 20 -P -Y -C /scratch/general/nfs1/u6071015/odgi -o ${SCAFFOLD}_sortedPGSGD.og

#odgi layout -i ${SCAFFOLD}_sorted.og -o ${SCAFFOLD}_sorted.lay -P --threads 20 

odgi viz \
    -i ${SCAFFOLD}_sortedPGSGD.og \
    -o ${SCAFFOLD}_sortedPGSGD.png \
    -x 4000 \
    -y 1500 
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
#a few modules missing fomr the sync
pip install networkx
pip install numpy
pip install uv
pip instal psutil
```
Okay, so now I am going to move to the /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree folder and write the python script from github and and sbatch script to run it. Here is the python script, which I will call pantree_config.py

```
from graph_var import PangenomeGraph

# Read a .gfa file
gfa_path = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/cactus/chrom-alignments/Scaffold_12__1_contigs__length_47609450.gfa "
#reference_path_index = 1
G, walks, walk_sample_names = PangenomeGraph.from_gfa_line_by_line(gfa_path, 
                                                return_walks=True, ref_name="Hap2_t_crist_hwy154_cen4119.2")
                                                #reference_path_index=reference_path_index)

# Generate vcf file
vcf_path = "/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/HWY154_REF_4119Hap2_pantree.vcf"
chr_id = "Scaffold_12"
G.write_vcf(gfa_path, vcf_path, chr_id)

# Enumerate variants of different types
edge_type_count: dict = G.variant_edges_summary()

# Get the genotype of a walk, then reconstruct edge visit counts
genotype: dict = G.genotype(walks[0])
edge_visit_counts: dict = G.count_edge_visits(genotype)

```
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
#SBATCH -e pantree-%j.err
#SBATCH -o pantree-%j.out

module load miniforge3
cd /uufs/chpc.utah.edu/common/home/u6071015/software/pantree/graph_var
source .venv/bin/activate
export PYTHONPATH="/uufs/chpc.utah.edu/common/home/u6071015/software/pantree:${PYTHONPATH}"

# Run pantree
python /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/pantree/pantree_config.py

```

## Pairwise comparison in Progressive Cactus

We are also going to call SVs from the pairwise comparisons, specifically focusing on comparisons between the reference haplotype used for the pangenome (HWY154 Stripe Haplotype2) and the other haplotypes. For this, I am creating softlinks in ''/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus'' to existing hal files, and need to make the hal file for H154 Stripe 2/Refugio Stripe 1 pair. ti do so, I run the script run_cactus.sh in the directory. The script also requires input file cactusTcrGSH2_TcrGSR1.txt

```
(TcrGSH2:0.010,TcrGSR1:0.010);

TcrGSH2 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_hwy154_cen4119_hap2.fasta.masked
TcrGSR1 /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/genomes/t_crist_refug_cen4122_hap1.fasta.masked
```

```
#!/bin/sh 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus
#SBATCH --qos gompert-grn
#SBATCH -e cactus-%j.err
#SBATCH -o cactus-%j.out

cd /scratch/general/nfs1/u6071015/cactusNp

module load cactus/1.0.0

cactus timemajobStore /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusTcrGSH2_TcrGSR1.txt cactusStripe_TcrGSH2_TcrGSR1.hal --maxCores 80

cp /scratch/general/nfs1/cactusNp/cactusStripe_TcrGSH2_TcrGSR1.hal /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus/cactusStripe_TcrGSH2_TcrGSR1.hal

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/progressive_cactus

#Summarize Mutations
halSummarizeMutations *.hal
```

