source("./functions/disease2frame.R")

require(readr)
require(dplyr)
require(biomaRt)
require(msigdbr)
require(igraph)

#Reading in biogrid data

links = read_tsv("./data/BIOGRID-ALL-3.5.174.mitab.txt")

links = links %>% filter(`Taxid Interactor A` == "taxid:9606" & `Taxid Interactor B` == "taxid:9606")

#Split based on Biogrid ID system then pull out gene names

protein1 = lapply(strsplit(links$`Alt IDs Interactor A`, "locuslink:"), function(x) strsplit(x[2], "\\|")[[1]][1])
protein2 = lapply(strsplit(links$`Alt IDs Interactor B`, "locuslink:"), function(x) strsplit(x[2], "\\|")[[1]][1])

protein1 = unlist(protein1)
protein2 = unlist(protein2)

#Place into matrix in form that disease2frame works on

links = links %>% mutate(protein1 =protein1, 
                         protein2 = protein2) %>% dplyr::select(protein1, protein2)

links = links %>% dplyr::select(hgnc_symbol_a = protein1,hgnc_symbol_b = protein2)

#Getting msigs 

msig = msigdbr()

msig = msig %>% filter(gs_subcat %in% c("CC", "CP:KEGG", "CP:REACTOME") )

mod_over = bind_rows()
for (j in unique(msig$gs_name)) {
  curr_msig = msig %>% filter(gs_name == j)
  write.table(curr_msig$gene_symbol, file = "temp_rea_bio.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  diseases_curr = diseases2frame("temp_rea_bio.txt", links = links)
  mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = diseases_curr)
  mod_slope = coef(mod)[[1]]
  mod_sum = data.frame(Disease = j, Count =nrow(diseases_curr), Slope = mod_slope,
                       Category = curr_msig$gs_subcat[1])
  mod_over = bind_rows(mod_over, mod_sum)
  saveRDS(mod_over,"KEGG_Biogrid.rds")
}

#Redoing with every category at once

msig = msigdbr()

mod_over = bind_rows()
for (j in unique(msig$gs_name)) {
  curr_msig = msig %>% filter(gs_name == j)
  write.table(curr_msig$gene_symbol, file = "temp_rea_bio.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  diseases_curr = diseases2frame("temp_rea_bio.txt", links = links)
  mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = diseases_curr)
  mod_slope = coef(mod)[[1]]
  mod_sum = data.frame(Disease = j, Count =nrow(diseases_curr), Slope = mod_slope,
                       Category = curr_msig$gs_subcat[1])
  print(mod_sum)
  mod_over = bind_rows(mod_over, mod_sum)
  saveRDS(mod_over,"All_Msig_Biogrid.rds")
}
