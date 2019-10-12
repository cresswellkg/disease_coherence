source("perm_slopes.R")
require(dplyr)
require(readr)
require(igraph)

network = read.csv("coherence_string_filt.csv")
Count_Frame = read.csv("string_filt_edges_new.csv")

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

string = string %>% filter(combined_score>=500)


Count_Frame = Count_Frame %>% dplyr::filter(Internal!=0) %>%
group_by(Disease) %>% mutate(Count = n()) %>% filter(Count>9)

over_p = rbind()
for (i in unique(Count_Frame$Disease)[31:40]) {
  print(i)
  perm = slope_perm(Count_Frame, Disease_Name = i, Norm_Reg =NULL, database = string, permutations = 10000, Disease = TRUE)
  curr_p = cbind(i, perm)
  over_p = rbind(over_p, curr_p)
  saveRDS(over_p, "string_filt_p_31_40.rds")
}

