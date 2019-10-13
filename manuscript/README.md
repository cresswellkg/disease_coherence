# README

- `manuscript.*` - RMarkdown version of the manuscript

- `Figures.Rmd` - RMarkdown for making a PPTX out of figures in `Figures` folder

## `Figures`

- `Fig1_network_distributions.png` - **Figure 1. Degree distributions are affected by network size and configuration.** Degree distributions of A) ring networks and B) fully connected networks containing 10 and 20 nodes. While degree distributions of similarly coherent ring networks can be directly compared, degree distributions of fully connected networks depend on network size.

- `Fig2_Size_Slope_Plot.png` - **Figure 2. Network size is inversely associated with coherence.** Slopes of internal vs. external degree distributions (aka coherence) for KEGG and random networks are plotted against network size.

- `Fig3_summary_traits.png` - **Figure 3. Coherence estimates of traits.** Size of dots represents the level of normalized coherence (X-axis) for individual traits (Y-axis).

- `Fig4_summary_diseases.png` - **Figure 4. Coherence estimates of diseases.** Size of dots represents the level of normalized coherence (X-axis) for individual diseases (Y-axis).

- `FigS1_corr_coherence.png` - **Supplementary Figure S1. Normalization of coherence estimates alleviates its network size dependence.** Red/blue gradient and numbers represent Pearson correlation coefficients for each pairwise comparison of the number of SNPs, genes, and the normalized coherence estimates using the corresponding PPI databases.

- `FigS2_summary_traits_Biogrid.png` - **Supplementary Figure S2. Coherence estimates of traits using Biogrid as a reference PPI database.** Size of dots represents the level of normalized coherence (X-axis) for individual traits (Y-axis). Plots are faceted by categories. Missing entries indicate that, for a given trait, a network could not be built and the coherence cannot be estimated.

- `FigS3_summary_traits_STRING.png` - **Supplementary Figure S3. Coherence estimates of traits using STRING as a reference of PPI database.** See legend for Supplementary Figure S2.

- `FigS4_summary_diseases_Biogrid.png` - **Supplementary Figure S4. Coherence estimates of diseases using Biogrid as a reference PPI database.** See legend for Supplementary Figure S2.

- `FigS5_summary_diseases_STRING.png` - **Supplementary Figure S5. Coherence estimates of diseases using STRING as a reference PPI database.** See legend for Supplementary Figure S2.


## `Tables` - Created by Kellen

- `table_1.xlsx` - **Table 1. Summary statistics of the analyzed categories.** Average coherence estimates are shown. "NA" values indicate that a category lacked phenotypes with a sufficient number of PPIs. Sorted by "String Filtered" average coherence.

- `supplementary_table_S1.xlsx` - **Supplementary Table S1. Summary statistics of the analyzed phenotype categories.** Minimum, mean, median, maximum of coherence estimates, the number of genes, SNPs, the total number of phenotype networks. Each worksheet corresponds to the PPI database used.

- `supplementary_table_S2.xlsx` - **Supplementary Table S2. Summary of networks from MSigDB categories.** Minimum, mean, median, maximum of the slopes, the number of genes, and the total number of networks are shown for each collection, along with the correlation between network size and slope. "CP:/KEGG/Reactome/Biocarta" - canonical pathways, "GOBP/GOMF/GOCC" - gene ontology biological processes, molecular functions, cellular component, "TFT" - transcription factor targets, "CGN" - cancer gene neighborhoods, "CM" - cancer modules, "MIR" - microRNA targets, "CGP" - chemical and genetic perturbations.

- `supplementary_table_S3.xlsx` - **Supplementary Table S3. Phenotype-specific coherence estimates.** The PPI database (STRING, STRING filtered, Biogrid) are specified in the corresponding column names. For each phenotype category, the results are sorted by the "String Filt Coherence" column in descending order. "NA" indicates the corresponding value cannot be estimated due to lack of sufficient number of genes annotated with PPIs.