# Pangenome 

## Pangenome construction
- ```run_cactus-pangenome.sh``` make the pangenome, requires ```HWY154_REF.txt``` file
- ```run_gretl_graphstats.sh``` to evaluate the pangenome

## Pangenome visualization
- ```odgi_draw.sh```
- ```run_odgi_viz_multiplescaff.sh```

## Pangenome SV calling
- ```run_pantree.sh``` runs pantree
- ```run_pantree_summary.sh``` modified from pantree manuscript, outputs summaries of all SV
  - ```pantree_summary.py```
- ```run_summarizepantreecoords.sh``` input is vcf subset to inversions from pantree, output is a .tsv file of all inversion in reference coordinate space
- ```run_invpg_summarizelengths.sh``` summarizes output from INVPG_annot
