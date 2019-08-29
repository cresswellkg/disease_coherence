source("./functions/disease2frame.R")

library(readr) #install.packages("readr")
library(dplyr) #install.packages("dplyr")
library(biomaRt) #BiocManager::install("biomaRt")
library(msigdbr) #BiocManager::install("msigdbr")
library(igraph) #install.packages("igraph")

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

# msig = msigdbr()
# 
# msig = msig %>% filter(gs_subcat %in% c("CC", "CP:KEGG", "CP:REACTOME") )
# 
# mod_over = bind_rows()
# for (j in unique(msig$gs_name)) {
#   curr_msig = msig %>% filter(gs_name == j)
#   write.table(curr_msig$gene_symbol, file = "temp_rea_bio.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
#   diseases_curr = diseases2frame("temp_rea_bio.txt", links = links)
#   mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = diseases_curr)
#   mod_slope = coef(mod)[[1]]
#   mod_sum = data.frame(Disease = j, Count =nrow(diseases_curr), Slope = mod_slope,
#                        Category = curr_msig$gs_subcat[1])
#   mod_over = bind_rows(mod_over, mod_sum)
#   saveRDS(mod_over,"KEGG_Biogrid.rds")
# }

disease_overall_biogrid     <- readRDS("biogrid_edges_new.rds")

#BIOGRID edge counts

dis_slopes_biogrid_count = disease_overall_biogrid %>% filter(Internal != 0) %>% filter(Category != "Random") %>% filter(Category != "KEGG") %>%
  group_by(Disease, Category) %>% mutate(Count = n()) %>% group_by(Disease, Category, Count) %>%
  do({
    mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = .)
    data.frame(Slope= coef(mod)[1])
  }) %>% arrange(Slope) %>% filter(Count>4)

#Calculating normalized coherence

# KEGG
msig = readRDS("KEGG_Biogrid.rds")
msig =msig %>% filter(Category == "CP:KEGG")
coherence_tab = bind_rows()
for (i in unique(dis_slopes_biogrid_count$Disease)) {
  print(i)
  dis = dis_slopes_biogrid_count %>% filter(Disease == i)
  msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
  msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]

  diseases = bind_rows()
  j = 1
    while(j<11) {
    diseases_curr = GeneSample(dis$Count, links, directed = FALSE)
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
      print(coherence_tab)
      saveRDS(coherence_tab, "./data/Coherence_Results/coherence_biogrid.rds")
}

# GOCC
msig = readRDS("KEGG_Biogrid.rds")
msig =msig %>% filter(Category == "CC") %>% filter(!is.na(Slope))
coherence_tab = bind_rows()
for (i in unique(dis_slopes_biogrid_count$Disease)) {
  print(i)
  dis = dis_slopes_biogrid_count %>% filter(Disease == i)
  msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
  msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]

  diseases = bind_rows()
  j = 1
    while(j<11) {
    diseases_curr = GeneSample(dis$Count, links, directed = FALSE)
    diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
    if (any(diseases_curr$Internal !=0)) {
      j = j+1
      diseases = bind_rows(diseases, diseases_curr)

    }
    }

      mod_rand = lm(External ~ Internal-1, data = diseases)
      rand_slope = coef(mod_rand)[[1]]

      coherence = (dis$Slope - rand_slope)/(median(msigs$Slope, na.rm = TRUE)-rand_slope)
      coherence = data.frame(Disease = dis$Disease, Coherence = coherence)

      coherence_tab = bind_rows(coherence, coherence_tab)
      saveRDS(coherence_tab, "./data/Coherence_Results/coherence_biogrid_CC.rds")
}


#Reactome
msig = readRDS("KEGG_Biogrid.rds")
msig =msig %>% filter(Category == "CP:REACTOME") %>% filter(!is.na(Slope))
coherence_tab = bind_rows()
for (i in unique(dis_slopes_biogrid_count$Disease)) {
  print(i)
  dis = dis_slopes_biogrid_count %>% filter(Disease == i)
  msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
  msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]

  diseases = bind_rows()
  j = 1
    while(j<11) {
    diseases_curr = GeneSample(dis$Count, links, directed = FALSE)
    diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
    if (any(diseases_curr$Internal !=0)) {
      j = j+1
      diseases = bind_rows(diseases, diseases_curr)

    }
    }

      mod_rand = lm(External ~ Internal-1, data = diseases)
      rand_slope = coef(mod_rand)[[1]]

      coherence = (dis$Slope - rand_slope)/(median(msigs$Slope, na.rm = TRUE)-rand_slope)
      coherence = data.frame(Disease = dis$Disease, Coherence = coherence)

      coherence_tab = bind_rows(coherence, coherence_tab)
      saveRDS(coherence_tab, "./data/Coherence_Results/coherence_biogrid_react.rds")
}


