# ToDo

- Figure 1 and 2

- Table 1 and all supplementary tables

- See `GWAS_Prep.Rmd` and `Connectivity_Comparison.Rmd` - organize your code like this. Key points - make more accurate indents, avoid extra lines, comment liberally, assume user has nothing installed. It saves time for both of us as I'll bother you less with such requests.

- `Coherence_Calculations_New.Rmd` - split commented code into a separate Rmd, make it run uncommented, make code release-ready, describe in README time it takes to run and the fact of pre-generated rds files

- Make other code release-ready, including functions

## Questions

- Is this description accurate?

- `Connectivity_Comparison.Rmd` - Calculation of  degree distributions and unnormalized coherence for disease related data. Includes old analyisis of random networks and KEGG pathways. Results are produced for Biogrid and STRING.

My understanding, it just makes PPI rds objects, no?


- Do we need, where it is, this file? - `efo.owl` - Experimental factor ontology file. Used to generate categories
    - What EFO code, currently commented, does? Current comments are unclear.

- In `Connectivity_Comparison.Rmd`, do we need all the functions loaded?

- How `KEGG_String.rds` and the like in "Analysis" folder were created?

