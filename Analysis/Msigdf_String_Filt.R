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

string_filt = string %>% filter(combined_score>=500)

disease_overall_string_filt <- readRDS("string_filt_edges_new.rds")

#STRING Filtered
dis_slopes_string_filt_count = disease_overall_string_filt %>% filter(Internal != 0) %>% filter(Category != "Random") %>% filter(Category != "KEGG") %>%
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
  print(j)
  curr_msig = msig %>% filter(gs_name == j)
  write.table(curr_msig$gene_symbol, file = "temp_rea_string_filt.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  diseases_curr = diseases2frame("temp_rea_string_filt.txt", links = string_filt)
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

saveRDS(mod_over, "KEGG_String_Filt.rds")

# KEGG
#Get msig pathway genes for KEGG, Reactome and GOCC
msig = readRDS("KEGG_String_Filt.rds")
#Get pathway specific counts
Count_Sum = msig %>% dplyr::select(Disease,Category, Num_Genes) %>% distinct()
#Splitting into list by categories
Count_Sum = split(Count_Sum, Count_Sum$Category)
coherence_tab = bind_rows()
for (i in unique(dis_slopes_string_filt_count$Disease)) {
  #Get disease specific genes
  dis = dis_slopes_string_filt_count %>% filter(Disease == i)
  #Get 30 closest msigdf networks taking 10 from each database and bind
  close_nets = bind_rows(lapply(Count_Sum,  function(x) x %>% mutate(Diff = abs(Num_Genes-dis$Count)) %>%
                                  arrange(Diff) %>% head(.,10)))
  #Subsetting msig to only include close pathways
  msig_sub = msig %>% filter(Disease %in% close_nets$Disease)
  #Square rooting because its not already square rooted
  mod_msig = lm(I(sqrt(External)) ~ Internal-1, data = msig_sub)
  #Pull out coefficient
  mod_slope = coef(mod_msig)[[1]]
  diseases = bind_rows()
  j = 1
  while(j<100) {
    diseases_curr = GeneSample(dis$Count, string_filt, directed = FALSE)
    diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
    if (any(diseases_curr$Internal !=0)) {
      j = j+1
      diseases = bind_rows(diseases, diseases_curr)
      
    }
  }
  
  mod_rand = lm(External ~ Internal-1, data = diseases)
  rand_slope = coef(mod_rand)[[1]]
  
  coherence = (dis$Slope - rand_slope)/(mod_slope-rand_slope)
  coherence = data.frame(Disease = dis$Disease, Coherence = coherence, 
                         Msig_Norm = mod_slope, Rand_Norm = rand_slope, Count = dis$Count)
  coherence_tab = bind_rows(coherence, coherence_tab)
  print(coherence_tab)
  saveRDS(coherence_tab, "./data/Coherence_Results/coherence_string_filt.rds")
}
