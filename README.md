# Timema-SV-Methods
Identifying how different methods of structural variant (SV) detection identify SVs, and whether Does SV size, frequency, or age have an effect on SV detection in different methods? Using data from this Gompert et al. 2025 (https://www.science.org/doi/10.1126/science.adp3745), and specifically focused on inversions and translocations.

## Pangenome Creation
Starting off with making pangenomes with 1) 4 hwy154 genomes and 2) 8 genomes (4 hwy154 and 4 refugio) in cactus minigraph.
|  POP  | STRIPE |   ID  |
|-------|--------|-------|
|HWY154 |striped |cen4119|
|HWY154 |  green |cen4280|
|REFUGIO|striped |cen4122|
|REFUGIO|  green |cen4120|

Working directory can be found in ```/uufs/chpc.utah.edu/common/home/gompert-group3/projects/timema_SVmethods/```. I copied the repeat masked genomes to subdirectory "genomes" in this directory.

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
