source("./functions/disease2frame.R")
source("./functions/GeneSample.R")

library(readr) #install.packages("readr")
library(dplyr) #install.packages("dplyr")
library(biomaRt) #BiocManager::install("biomaRt")
library(msigdbr) #BiocManager::install("msigdbr")
library(igraph) #install.packages("igraph")

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

#Getting msigs 

# msig = msigdbr()
# 
# msig = msig %>% filter(gs_subcat %in% c("CC", "CP:KEGG", "CP:REACTOME") )
# 
# mod_over = bind_rows()
# for (j in unique(msig$gs_name)) {
#   curr_msig = msig %>% filter(gs_name == j)
#   write.table(curr_msig$gene_symbol, file = "temp_rea_string.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
#   diseases_curr = diseases2frame("temp_rea_string.txt", links = string)
#   if (nrow(diseases_curr) == 0) {
#     next
#   }
#   mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = diseases_curr)
#   mod_slope = coef(mod)[[1]]
#   mod_sum = data.frame(Disease = j, Count =nrow(diseases_curr), Slope = mod_slope,
#                        Category = curr_msig$gs_subcat[1])
#   mod_over = bind_rows(mod_over, mod_sum)
#   saveRDS(mod_over,"KEGG_String.rds")
# }

#Redoing with every category at once

msig = msigdbr()

msig = msig %>% filter(gs_subcat != "")
mod_over = bind_rows()
for (j in unique(msig$gs_name)) {
  curr_msig = msig %>% filter(gs_name == j)
  write.table(curr_msig$gene_symbol, file = "temp_rea_string.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  diseases_curr = diseases2frame("temp_rea_string.txt", links = string)
  if (nrow(diseases_curr) == 0) {
    next
  }
  mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = diseases_curr)
  mod_slope = coef(mod)[[1]]
  mod_sum = data.frame(Disease = j, Count =nrow(diseases_curr), Slope = mod_slope,
                       Category = curr_msig$gs_subcat[1])
  mod_over = bind_rows(mod_over, mod_sum)
  saveRDS(mod_over,"./data/Coherence_Results/All_Msig_String.rds")
}

#Calculating normalized coherence

disease_overall_string  <- readRDS("string_edges_new.rds")

#Reading in slopes for diseases

#Calculate slopes for all diseases under each category
dis_slopes_string_count = disease_overall_string %>% filter(Internal != 0) %>% filter(Category != "Random") %>% filter(Category != "KEGG") %>%
  group_by(Disease, Category) %>% mutate(Count = n()) %>% group_by(Disease, Category, Count) %>%
  do({
    mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = .)
    data.frame(Slope= coef(mod)[1])
  }) %>% arrange(Slope) %>% filter(Count>4)

#Read in KEGG pathways for coherence
msig = readRDS("KEGG_String.rds")
msig =msig %>% filter(Category == "CP:KEGG")
coherence_tab = bind_rows()
for (i in unique(dis_slopes_string_count$Disease)) {
  #Filter to only include specific disease coherence
  dis = dis_slopes_string_count %>% ungroup() %>% filter(Disease == i)
  #Get closest KEGG pathways in term of gene counts
  msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
  #Take top 10 KEGG pathways
  msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]

  diseases = bind_rows()
  #Generate random networks
  j = 1
    while(j<11) {
    #Sample gene networks equal in size to disease count
    diseases_curr = GeneSample(dis$Count, string, directed = FALSE)
    #Get random gene networks
    diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))

    #Repeat if there are no internal edges
    if (any(diseases_curr$Internal !=0)) {
      j = j+1
      diseases = bind_rows(diseases, diseases_curr)

    }
    }
    #Calculate overall slope of random disease networks
      mod_rand = lm(External ~ Internal-1, data = diseases)
      rand_slope = coef(mod_rand)[[1]]

      #Normalize based on random slope and the median of kegg slopes (More robust)
      coherence = (dis$Slope - rand_slope)/(median(msigs$Slope, na.rm = TRUE)-rand_slope)
      coherence = data.frame(Disease = dis$Disease, Coherence = coherence, Msig_Norm = median(msigs$Slope, na.rm = TRUE), Rand_Norm = rand_slope )
      coherence_tab = bind_rows(coherence, coherence_tab)
      saveRDS(coherence_tab, "./data/Coherence_Results/coherence_string.rds")
}

#Reactome

msig = readRDS("KEGG_String.rds")
msig =msig %>% filter(Category == "CP:REACTOME") %>% filter(!is.na(Slope))
coherence_tab = bind_rows()
for (i in unique(dis_slopes_string_count$Disease)) {
  print(i)
  dis = dis_slopes_string_count %>% filter(Disease == i)
  msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
  msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]

  diseases = bind_rows()
  j = 1
    while(j<11) {
    diseases_curr = GeneSample(dis$Count, string, directed = FALSE)
    diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
    if (any(diseases_curr$Internal !=0)) {
      j = j+1
      diseases = bind_rows(diseases, diseases_curr)

    }
    }

      mod_rand = lm(External ~ Internal-1, data = diseases)
      rand_slope = coef(mod_rand)[[1]]

      coherence = (dis$Slope - rand_slope)/(median(msigs$Slope, na.rm = TRUE)-rand_slope)
      coherence = data.frame(Disease = dis$Disease, Coherence = coherence, Msig_Norm = median(msigs$Slope, na.rm = TRUE), Rand_Norm = rand_slope )

      coherence_tab = bind_rows(coherence, coherence_tab)
      saveRDS(coherence_tab, "./data/Coherence_Results/coherence_string_react.rds")
}



# CC

msig = readRDS("KEGG_String.rds")
msig =msig %>% filter(Category == "CC") %>% filter(!is.na(Slope))
coherence_tab = bind_rows()
for (i in unique(dis_slopes_string_count$Disease)) {
  print(i)
  dis = dis_slopes_string_count %>% filter(Disease == i)
  msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
  msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]

  diseases = bind_rows()
  j = 1
    while(j<11) {
    diseases_curr = GeneSample(dis$Count, string, directed = FALSE)
    diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
    if (any(diseases_curr$Internal !=0)) {
      j = j+1
      diseases = bind_rows(diseases, diseases_curr)

    }
    }

      mod_rand = lm(External ~ Internal-1, data = diseases)
      rand_slope = coef(mod_rand)[[1]]

      coherence = (dis$Slope - rand_slope)/(median(msigs$Slope, na.rm = TRUE)-rand_slope)
      coherence = data.frame(Disease = dis$Disease, Coherence = coherence, Msig_Norm = median(msigs$Slope, na.rm = TRUE), Rand_Norm = rand_slope )

      coherence_tab = bind_rows(coherence, coherence_tab)
      saveRDS(coherence_tab, "./data/Coherence_Results/coherence_string_CC.rds")
}

# Random networks

#Read in KEGG pathways for coherence
msig = readRDS("KEGG_String.rds")
msig = msig %>% filter(Category == "CP:KEGG")
string_random = bind_rows()
for (i in msig$Count) {
  j = 0
  diseases = bind_rows()
    while(j<11) {
    diseases_curr = GeneSample(i, string, directed = FALSE)
    diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
    if (any(diseases_curr$Internal !=0)) {
      j = j+1
      diseases = bind_rows(diseases, diseases_curr)

    }
    }

      mod_rand = lm(External ~ Internal-1, data = diseases)
      diseases_curr = data.frame(Count = i, Slope = coef(mod_rand)[1])
      string_random = bind_rows(string_random, diseases_curr)
      print(string_random)
      saveRDS(string_random, "./data/Coherence_Results/string_random_edges.rds")
}


