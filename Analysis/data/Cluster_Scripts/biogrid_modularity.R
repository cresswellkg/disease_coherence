source("perm_slopes.R")
require(dplyr)
require(readr)
require(igraph)
source("./functions/diseases2mod.R")

network = read.csv("coherence_biogrid.csv")
Count_Frame = read.csv("biogrid_edges_new.csv")

Cat_Frame = Count_Frame %>% dplyr::select(Disease, Category) %>% distinct()

network = left_join(network, Cat_Frame)


links = read_tsv("BIOGRID-ALL-3.5.174.mitab.txt")

links = links %>% filter(`Taxid Interactor A` == "taxid:9606" & `Taxid Interactor B` == "taxid:9606")

#Split based on Biogrid ID system then pull out gene names

protein1 = lapply(strsplit(links$`Alt IDs Interactor A`, "locuslink:"), function(x) strsplit(x[2], "\\|")[[1]][1])
protein2 = lapply(strsplit(links$`Alt IDs Interactor B`, "locuslink:"), function(x) strsplit(x[2], "\\|")[[1]][1])

protein1 = unlist(protein1)
protein2 = unlist(protein2)

#Place into matrix in form that disease2frame works on

links = links %>% mutate(protein1 =protein1, 
                         protein2 = protein2) %>% dplyr::select(protein1, protein2)

links = links %>% rename(hgnc_symbol_a = protein1) %>% rename(hgnc_symbol_b = protein2)

Count_Frame = Count_Frame %>% dplyr::filter(Internal!=0) %>%
group_by(Disease) %>% mutate(Count = n()) %>% filter(Count>9)

over_mod = rbind()
for (i in unique(Count_Frame$Disease)) {
  Count_Curr = Count_Frame %>% filter(i %in% Disease)
   Curr_Net = links %>% filter( (hgnc_symbol_a %in% Count_Curr$Genes) | (hgnc_symbol_b %in% Count_Curr$Genes))
  mod = diseases2mod(Count_Curr, links)
  curr_mod = cbind(i, mod)
  over_mod = rbind(over_mod, curr_mod)
  saveRDS(over_mod, "biogrid_mod.rds")
}