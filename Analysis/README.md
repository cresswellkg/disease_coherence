# README

- `GWAS_Prep.Rmd` - Code for parsing through the ebicat and KEGG databases and creating gene lists for each disease and pathway. 
- `Connectivity_Comparison_Undirected.Rmd` - Calculation of  degree distributions and unnormalized coherence for disease related data. Includes old analyisis of random networks and KEGG pathways. Results are produced for Biogrid and STRING without filtering. Relies on `GWAS_Prep.Rmd` 
- `Connectivity_Comparison_Filtered.Rmd` - Similar to comparison but calculates unnormalized coherence for STRING but with scores below 500 filtered. Includes old analyisis of random networks and KEGG pathways. Relies on `GWAS_Prep.Rmd` 
-`Coherence_Calculations_New.Rmd` - Calculates normalized coherence using Msigdf and random networks. Relies on Msigdf.rds, Connectivity_Comparison_Filtered.Rmd and Connectivity_Comparison_Undirected.Rmd

# Folders

## functions
Contains functions for conversion of gene networks to degree distributions

## Data
Contains the data necessary to run the scripts in the folder

- `get_data.sh` - Shell script for downloading Biogrid and STRING PPI data
- `efo.owl` - Experimental factor ontology file. Used to generate categories
- 9606.protein.links.v11.0.txt.gz - STRING PPI database data. Needed for Connectivity_Comparison_*. Produced from get_data.sh
- BIOGRID-ALL-3.5.174.mitab.txt - Biogrid PPI database data. Needed for Connectivity_Comparison_undirected. Produced from running get_data.sh and unzipping BIOGRID-ALL-3.5.174.mitab.zip




