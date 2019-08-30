# README

## Instructions to reproduce the results

- Download the data 
````
cd Analysis/data
chmod +x get_data.sh
./get_data.sh
cd ..
````

- In RStudio, set working directory to the `Analysis` folder (`setwd("Analysis")`) and run the scripts in the following order:
    - `GWAS_Prep.Rmd` - Code for parsing through the ebicat and KEGG databases and creating gene lists for each disease and pathway. Creates `data/Disease_Genes` folder. 
    - `Connectivity_Comparison.Rmd` - Calculation of degree distributions for disease related data. Results are produced for Biogrid and STRING. Creates `biogrid_edges_new.rds`, `string_edges_new.rds`, `string_filt_edges_new.rds`.
    -`Coherence_Calculations_New.Rmd` - Takes normalized coherence, p-values and SNP/Gene counts and puts them into tables. Relies on coherence_*.rds files from Msigdf_*.R, `Connectivity_Comparison.Rmd` and supplementary_table_diseases_selected.csv. Outputs "./manuscript/tables/supplementary_table_S3.csv".
    -`Category_Tables.Rmd` - Produces tables and plots summarizing coherence of Msigdf pathways and diseases seperated by coherence. File produces "./manuscript/tables/Table_1.csv" and "./manuscript/tables/supplementary_table_S1.xlsx" and "Size_Slope_Plot.png". Relies on "./manuscript/tables/supplementary_table_S3.csv" from `Coherence_Calculations_New.Rmd`.

- `Msigdf_*.R` - Scripts for calculating normalized coherence for KEGG, REACTOME and CC networks. Each one outputs a dataset containing slopes, counts, pathway name and category for STRING, STRING Filtered and Biogrid and three normalized coherence files for KEGG, Reactome and GOCC. `Msigdf_String.R` produces two additional files `All_Msig_String.rds` containing slopes for all Msigdf pathways and `string_random_edges.rds` containing slopes for randomly generated pathways. Produce coherence_*.rds files that are used in Coherence_Calculations_New.Rmd

Scripts must be ran in the following order:

`GWAS_Prep.Rmd` -> `Connectivity_Comparison.Rmd` -> `Msigdf_*.R` (Any order) -> `Coherence_Calculations_New.Rmd` -> `Category_Tables.Rmd`

# Folders

## functions
Contains functions for conversion of gene networks to degree distributions

## Data
Contains the data necessary to run the scripts in the folder

- `get_data.sh` - Shell script for downloading Biogrid and STRING PPI data
- `efo.owl` - Experimental factor ontology file. Used to generate categories
- 9606.protein.links.v11.0.txt.gz - STRING PPI database data. Needed for Connectivity_Comparison_*. Produced from get_data.sh
- BIOGRID-ALL-3.5.174.mitab.txt - Biogrid PPI database data. Needed for Connectivity_Comparison_undirected. Produced from running get_data.sh and unzipping BIOGRID-ALL-3.5.174.mitab.zip

##Tables

- supplementary_table_S1.xlsx - Contains SNP/Gene counts, coherence and permutation p-values for STRING, BioGrid and STRING Filtered disease networks




