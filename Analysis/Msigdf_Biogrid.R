source("./functions/disease2frame.R")
source("./functions/GeneSample.R")

library(readr) #install.packages("readr")
library(dplyr) #install.packages("dplyr")
library(biomaRt) #BiocManager::install("biomaRt")
library(msigdbr) #BiocManager::install("msigdbr")
library(igraph) #install.packages("igraph")
library(tidyr) #install.packages("tidyr")
library(broom) #install.packages("broom")

#Option for number of msig pathways to take from each category
num_msig = 10
#Option for number of random networks to generate (Will generate num_rand-1 networks)
num_rand = 11


#Reading in biogrid data

links = read_tsv("./data/BIOGRID-ALL-3.5.174.mitab.txt.gz")

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

mod_over = list()
index = 1
for (j in unique(msig$gs_name)) {
  print(j)
  curr_msig = msig %>% filter(gs_name == j)
  write.table(curr_msig$gene_symbol, file = "temp_rea_bio.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  diseases_curr = diseases2frame("temp_rea_bio.txt", links = links)
  diseases_curr = diseases_curr %>% mutate(Disease = j)
  mod_over[[index]] = diseases_curr
  index =  index+1
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

saveRDS(mod_over, "KEGG_Biogrid.rds")

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
#Get msig pathway genes for KEGG, Reactome and GOCC
msig = readRDS("KEGG_Biogrid.rds")
#Get pathway specific counts
Count_Sum = msig %>% dplyr::filter(Internal>0) %>% 
  group_by(Disease) %>%
  mutate(Num_Genes = n()) %>% dplyr::select(Disease,Category, Num_Genes) %>% distinct()
#Splitting into list by categories
Count_Sum = split(Count_Sum, Count_Sum$Category)
coherence_tab = bind_rows()
for (i in unique(dis_slopes_biogrid_count$Disease)) {
  #Get disease specific genes
  dis = dis_slopes_biogrid_count %>% filter(Disease == i)
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
    diseases_curr = GeneSample(dis$Count, links, directed = FALSE)
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
      saveRDS(coherence_tab, "./data/Coherence_Results/coherence_biogrid.rds")
}



