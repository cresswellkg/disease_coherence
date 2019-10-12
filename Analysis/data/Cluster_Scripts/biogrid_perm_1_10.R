source("perm_slopes.R")
require(dplyr)
require(readr)
require(igraph)

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

over_p = rbind()
for (i in unique(Count_Frame$Disease)[1:10]) {
  print(i)
  perm = slope_perm(Count_Frame, i, Norm_Reg =NULL, database = links, permutations = 10000, Disease = TRUE)
  curr_p = cbind(i, perm)
  over_p = rbind(over_p, curr_p)
  saveRDS(over_p, "biogrid_p_1_10.rds")
}

