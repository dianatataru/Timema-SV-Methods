# Comparative Alignment

## Whole-genome alignment and visualization with Cactus
- ```run_cactus.sh```
This uses cactus version 2.7.2 to align whole-genome assemblies pairwise and create .hal files
- ```runHalSynteny.sh```
This takes the .hal files and makes .psl files for visualization
- ```plottingalignment.R```
This takes the .psl files and makes pairwise dotplots.

## SV calling with SyRI
You can use the above progressive cactus runs to figure out which chromosomes to invert and fuse. SyRI requires chromosomes to be in the same configuration.
- ```invert_chroms.py```
Invert chromosomes that are oriented differently in the alignments
- ```create_fused_references.py```
Fuse chromosomes in the reference in multiple versions to match those that are fused on the other genomes
- ```rename_fasta_headers_withfused.py```
To rename the chromosomes to match the reference, and give fused versions unique names
- ```run_syRI.sh```
Finally you can run syRI to call SV, pairwise, across the genomes. 
