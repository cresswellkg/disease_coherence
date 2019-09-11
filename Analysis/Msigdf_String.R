source("./functions/disease2frame.R")
source("./functions/GeneSample.R")

library(readr) #install.packages("readr")
library(dplyr) #install.packages("dplyr")
library(biomaRt) #BiocManager::install("biomaRt")
library(msigdbr) #BiocManager::install("msigdbr")
library(igraph) #install.packages("igraph")

#Option for number of msig pathways to take from each category
num_msig = 10
#Option for number of random networks to generate (Will generate num_rand-1 networks)
num_rand = 11

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

#Redoing with every category at once

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

#STRING Filtered
dis_slopes_string_count = disease_overall_string %>% filter(Internal != 0) %>% 
  filter(Category != "Random") %>% filter(Category != "KEGG") %>%
  group_by(Disease, Category) %>% mutate(Count = n()) %>% group_by(Disease, Category, Count) %>%
  do({
    mod = lm(sqrt(External) ~ sqrt(Internal)-1, data = .)
    data.frame(Slope= coef(mod)[1])
  }) %>% arrange(Slope) %>% filter(Count>4)

#Getting msigs 

msig = msigdbr()

msig = msig %>% filter(gs_subcat %in% c("CC", "CP:KEGG", "CP:REACTOME") )

mod_over = list()
index = 1
for (j in unique(msig$gs_name)) {
  curr_msig = msig %>% filter(gs_name == j)
  write.table(curr_msig$gene_symbol, file = "temp_rea_string.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  diseases_curr = diseases2frame("temp_rea_string.txt", links = string)
  diseases_curr = diseases_curr %>% mutate(Disease = j)
  mod_over[[index]] = diseases_curr
  index = index+1
}

mod_over = bind_rows(mod_over)

#Applying categories

msig_cats  = msig %>% 
  dplyr::select(Disease = gs_name, Category = gs_subcat) %>% 
  distinct() #Get categories in order

#Joining categories

mod_over = left_join(mod_over, msig_cats)

#Getting gene counts

mod_over = mod_over %>% group_by(Disease) %>% mutate(Num_Genes = n())

#Saving

saveRDS(mod_over, "KEGG_String.rds")

#Calculating normalized coherence

# KEGG
#Get msig pathway genes for KEGG, Reactome and GOCC
msig = readRDS("KEGG_String.rds")
#Get pathway specific counts
Count_Sum = msig %>% dplyr::filter(Internal>0) %>% 
  group_by(Disease) %>%
  mutate(Num_Genes = n()) %>% dplyr::select(Disease,Category, Num_Genes) %>% distinct()
#Splitting into list by categories
Count_Sum = split(Count_Sum, Count_Sum$Category)
coherence_tab = bind_rows()
for (i in unique(dis_slopes_string_count$Disease)) {
  #Get disease specific genes
  dis = dis_slopes_string_count %>% filter(Disease == i)
  #Get 30 closest msigdf networks taking 10 from each database and bind
  close_nets = bind_rows(lapply(Count_Sum,  function(x) x %>% mutate(Diff = abs(Num_Genes-dis$Count)) %>%
                                  arrange(Diff) %>% head(.,num_msig)))
  #Subsetting msig to only include close pathways
  msig_sub = msig %>% filter(Disease %in% close_nets$Disease)
  #Taking slopes for each category
  mod_msig = msig_sub %>% group_by(Disease) %>% do(tidy(
    coef(lm(I(sqrt(External)) ~ I(sqrt(Internal))-1, data = .))[[1]]
  ))
  #Take median
  mod_slope = median(mod_msig$x)
  #Creating 10 random networks
  j = 1
  slope_list = list()
  while(j<num_rand) {
    #Empty data frame for disease networks
    diseases = bind_rows()
    while(nrow(diseases)<dis$Count) {
      #Generate random diseases network
      diseases_curr = GeneSample(dis$Count, string, directed = FALSE)
      #Save internal and external edges
      diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = diseases_curr$internal_connectivity, External = diseases_curr$external_connectivity)
      #Filter out zeros
      diseases_curr = diseases_curr %>% filter(Internal !=0)
      #Bind together non-zero internal/external edges
      diseases = bind_rows(diseases, diseases_curr)
    }
    #Save slope of current network
    slope_list[[j]] = coef(lm(I(sqrt(External)) ~ I(sqrt(Internal))-1, data = diseases))[[1]]
    #Iterate
    j = j+1
    
  }
  
  rand_slope = median(unlist(slope_list))
  
  coherence = (dis$Slope - rand_slope)/(mod_slope-rand_slope)
  coherence = data.frame(Disease = dis$Disease, Coherence = coherence, 
                         Msig_Norm = mod_slope, Rand_Norm = rand_slope, Count = dis$Count)
  coherence_tab = bind_rows(coherence, coherence_tab)
  print(coherence_tab)
  saveRDS(coherence_tab, "./data/Coherence_Results/coherence_string.rds")
}



