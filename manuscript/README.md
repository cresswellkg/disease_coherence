# README

- `Analysis_Tables_New.Rmd` - analysis of new tables. New analysis refers to disease-matched random and KEGG references.

## `Figures`

- `Fig1_network_distributions.png` - Created by Kellen
- `Size_Slope_Plot.png` - Created by Kellen

- `Fig3_summary_traits.png` - Created by `Analysis_Tables_New.Rmd`
- `Fig4_summary_diseases.png` - Created by `Analysis_Tables_New.Rmd`
- `FigS1_corr_coherence.png` - Created by `Analysis_Tables_New.Rmd`
- `FigS2_corr_references.png` - Created by `Analysis_Tables_New.Rmd`
- `FigS3_summary_traits_KEGG.png` - Created by `Analysis_Tables_New.Rmd`
- `FigS4_summary_traits_GOCC.png` - Created by `Analysis_Tables_New.Rmd`
- `FigS5_summary_traits_Reactome.png` - Created by `Analysis_Tables_New.Rmd`
- `FigS6_summary_diseases_KEGG.png` - Created by `Analysis_Tables_New.Rmd`
- `FigS7_summary_diseases_GOCC.png` - Created by `Analysis_Tables_New.Rmd`
- `FigS8_summary_diseases_Reactome.png` - Created by `Analysis_Tables_New.Rmd`


## `Tables` - Created by Kellen

**Table 1.[FILE: supplementary_table_categories_short.csv]** Summary statistics of the analyzed categories. Average coherence estimates are shown. "NA" values indicate that a phenotype was not considered due to low number of genes with PPIs.

**Supplementary Table S1.[FILE: tables/supplementary_table_categories_new.xlsx] Summary statistics of the analyzed phenotypes.** Each worksheet corresponds to the PPI database used.

**Supplementary Table S2.[FILE: tables/Msigdf_Summary.xlsx] Summary of networks from MSigDB categories.** Minimum, maximum, mean and median of the slopes and number of genes are shown for each collection, along with the total number of networks and the correlation between network size and slope.

**Supplementary Table S3.[FILE: tables/supplementary_table_S1.csv] Phenotype-specific coherence estimates.** The PPI database (STRING, STRING filtered, Biogrid) and the reference of high coherence (KEGG, GOCC, Reactome) are specified in the corresponding columns for coherence and permutation p-value estimations.

- `tables/Mod_Table.xlsx` - same as `supplementary_table_S1.csv`, includes modularity



# Unneeded

- `DREAM_gene_extraction.Rmd` - Exploring gene summary statistics

## `data`

- `media-1.xlsx` - Trait-associated lists of genes were downloaded from the Disease Module Identification DREAM Challange (Choobdar et al. 2018). A comprehensive collection of 180 GWAS datasets, ranging over diverse disease-related human phenotypes, were used to annotate genes (Table 1). [Supplemental Table 1](https://www.biorxiv.org/content/10.1101/265553v2.supplementary-material)
    - Table 1. Sources for the 180 trait-specific networks. Columns: GWAS Set - indicator denoting whether the GWAS network was used for the leaderboard phase (Leaderboard) or the final evaluation (Final); Category - trait category; Trait-Type - disease trait; GWAS Type - GWAS specific disease trait; SNPs - total number of SNPs associated with the network; Study - name of study; Reference- reference of the study; Source - study URL; Data File - name of network dataset.


## `misc`

- `07-07-2019_Dozmorov_progress.Rmd` - brief progress report
- `manuscript 7-29_19 KSK.docx` - Mikhail, Here is the MS with comments. Overall, the writing was very good and clear. The major issue for me was what to say about the relationship between our measures of coherence and the underlying biology of the different traits. I think you go too far, especially in the results section making claims that are tenuous. I would favor just giving the results in the results section and then have a general paragraph in the discussion noting how tentative these finding are and uncertain our interpretations can be at this stage. Then you might give some broad trends – e.g. that on average more physiological complex traits/disorders tend to have lower coherence. But the relationship seems to me to be far from complete. I would try to avoiding getting into the details e.g. one kind of cancer versus another or trying to say the schizophrenia is more complex than bipolar illness or autism. Just too tenuous. Otherwise, quite an accomplishment. Ken. 2019-07-29

- `manuscript 7-29_19 KSKCFCcomments.docx` - I have some minor suggestions and a major one at the end. 2019-31-07
