source("perm_slopes.R")
require(dplyr)
require(readr)
require(igraph)
source("./functions/diseases2mod.R")

network = read.csv("coherence_string.csv")
Count_Frame = read.csv("string_edges_new.csv")

Cat_Frame = Count_Frame %>% dplyr::select(Disease, Category) %>% distinct()

network = left_join(network, Cat_Frame)

string = read_table2("9606.protein.links.v11.0.txt.gz")

string = string %>% mutate(protein1 = gsub("9606.", "", protein1), 
                           protein2 = gsub("9606.", "", protein2))

#Read in ensemble protein ID

genes = readRDS("genes.rds")

#Join gene names with string

string = left_join(string, genes, by = c("protein1" = "ensembl_peptide_id", "protein2" = "ensembl_peptide_id"))

string = string %>% rename(hgnc_symbol_a = hgnc_symbol)

string = left_join(string, genes, by = c("protein2" = "ensembl_peptide_id"))

string = string %>% rename(hgnc_symbol_b = hgnc_symbol)


Count_Frame = Count_Frame %>% dplyr::filter(Internal!=0) %>%
group_by(Disease) %>% mutate(Count = n()) %>% filter(Count>9)

over_mod = rbind()
for (i in unique(Count_Frame$Disease)) {
  Count_Curr = Count_Frame %>% filter(i %in% Disease)
   Curr_Net = string %>% filter( (hgnc_symbol_a %in% Count_Curr$Genes) | (hgnc_symbol_b %in% Count_Curr$Genes))
  mod = diseases2mod(Count_Curr, string)
  curr_mod = cbind(i, mod)
  over_mod = rbind(over_mod, curr_mod)
  saveRDS(over_mod, "string_mod.rds")
}
