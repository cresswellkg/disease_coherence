source("./functions/disease2frame.R")

require(readr)
require(dplyr)
require(biomaRt)
require(msigdbr)
require(igraph)

string = read_table2("./data/9606.protein.links.v11.0.txt.gz")

string = string %>% mutate(protein1 = gsub("9606.", "", protein1), 
                           protein2 = gsub("9606.", "", protein2))

#Read in ensemble protein ID

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "www")

genes = biomaRt::getBM(attributes=c('ensembl_peptide_id','hgnc_symbol'), mart = ensembl)


#Join gene names with string

string = left_join(string, genes, by = c("protein1" = "ensembl_peptide_id", "protein2" = "ensembl_peptide_id"))

string = left_join(string, genes, by = c("protein2" = "ensembl_peptide_id"))

string = string %>% dplyr::rename(hgnc_symbol_a = hgnc_symbol.x)

string = string %>% dplyr::rename(hgnc_symbol_b = hgnc_symbol.y)

string = string %>% filter(combined_score>=500)


#Getting msigs 

msig = msigdbr()

msig = msig %>% filter(gs_subcat %in% c("CC", "CP:KEGG", "CP:REACTOME") )

mod_over = bind_rows()
for (j in rev(unique(msig$gs_name))) {
  curr_msig = msig %>% filter(gs_name == j)
  write.table(curr_msig$gene_symbol, file = "temp_rea.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  diseases_curr = diseases2frame("temp_rea.txt", links = string)
  if (nrow(diseases_curr) == 0) {
    next
  }
  mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = diseases_curr)
  mod_slope = coef(mod)[[1]]
  mod_sum = data.frame(Disease = j, Count =nrow(diseases_curr), Slope = mod_slope,
                       Category = curr_msig$gs_subcat[1])
  print(mod_sum)
  mod_over = bind_rows(mod_over, mod_sum)
  saveRDS(mod_over,"KEGG_String_Filt.rds")
  }

#Redoing with every category at once

msig = msigdbr()

mod_over = bind_rows()
for (j in rev(unique(msig$gs_name))) {
  curr_msig = msig %>% filter(gs_name == j)
  write.table(curr_msig$gene_symbol, file = "temp_rea.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  diseases_curr = diseases2frame("temp_rea.txt", links = string)
  if (nrow(diseases_curr) == 0) {
    next
  }
  mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = diseases_curr)
  mod_slope = coef(mod)[[1]]
  mod_sum = data.frame(Disease = j, Count =nrow(diseases_curr), Slope = mod_slope,
                       Category = curr_msig$gs_subcat[1])
  print(mod_sum)
  mod_over = bind_rows(mod_over, mod_sum)
  saveRDS(mod_over,"All_Msig_String_Filt.rds")
}