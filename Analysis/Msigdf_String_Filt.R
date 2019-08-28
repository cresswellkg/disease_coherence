source("./functions/disease2frame.R")

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

# #Normalized coherence
# #KEGG
# msig = readRDS("KEGG_String_Filt.rds")
# msig =msig %>% filter(Category == "CP:KEGG")
# coherence_tab = bind_rows()
# for (i in unique(dis_slopes_string_filt_count$Disease)) {
#   print(i)
#   dis = dis_slopes_string_filt_count %>% filter(Disease == i)
#   msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
#   msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]
# 
#   diseases = bind_rows()
#   j = 1
#     while(j<11) {
#     diseases_curr = GeneSample(dis$Count, string_filt, directed = FALSE)
#     diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
#     if (any(diseases_curr$Internal !=0)) {
#       j = j+1
#       diseases = bind_rows(diseases, diseases_curr)
# 
#     }
#     }
# 
#       mod_rand = lm(External ~ Internal-1, data = diseases)
#       rand_slope = coef(mod_rand)[[1]]
# 
#       coherence = (dis$Slope - rand_slope)/(median(msigs$Slope, na.rm = TRUE)-rand_slope)
#     coherence = data.frame(Disease = dis$Disease, Coherence = coherence, Msig_Norm = median(msigs$Slope, na.rm = TRUE), Rand_Norm = rand_slope )
#       coherence_tab = bind_rows(coherence, coherence_tab)
#       saveRDS(coherence_tab, "./data/Coherence_Results/coherence_string_filt.rds")
# }

# #Reactome
# msig = readRDS("KEGG_String_Filt.rds")
# msig =msig %>% filter(Category == "CP:REACTOME") %>% filter(!is.na(Slope))
# coherence_tab = bind_rows()
# for (i in unique(dis_slopes_string_filt_count$Disease)) {
#   print(i)
#   dis = dis_slopes_string_filt_count %>% filter(Disease == i)
#   msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
#   msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]
#   
#   diseases = bind_rows()
#   j = 1
#     while(j<11) {
#     diseases_curr = GeneSample(dis$Count, string_filt, directed = FALSE)
#     diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
#     if (any(diseases_curr$Internal !=0)) {
#       j = j+1
#       diseases = bind_rows(diseases, diseases_curr)
# 
#     }
#     }
#   
#       mod_rand = lm(External ~ Internal-1, data = diseases)
#       rand_slope = coef(mod_rand)[[1]]
#       
#       coherence = (dis$Slope - rand_slope)/(median(msigs$Slope, na.rm = TRUE)-rand_slope)
#       coherence = data.frame(Disease = dis$Disease, Coherence = coherence)
#       
#       coherence_tab = bind_rows(coherence, coherence_tab)
#       saveRDS(coherence_tab, "./data/Coherence_Results/coherence_string_filt_react.rds")
# }
#   
#   
# 


# #GOCC
# msig = readRDS("KEGG_String_Filt.rds")
# msig =msig %>% filter(Category == "CC") %>% filter(!is.na(Slope))
# coherence_tab = bind_rows()
# for (i in unique(dis_slopes_string_filt_count$Disease)) {
#   print(i)
#   dis = dis_slopes_string_filt_count %>% filter(Disease == i)
#   msigs = msig %>% mutate(Num_Genes = abs(Count-dis$Count)) %>% arrange(Num_Genes)
#   msigs = msigs[msigs$Num_Genes %in% unique(msigs$Num_Genes)[1:10],]
# 
#   diseases = bind_rows()
#   j = 1
#     while(j<11) {
#     diseases_curr = GeneSample(dis$Count, string_filt, directed = FALSE)
#     diseases_curr = data.frame(Gene = names(diseases_curr$internal_connectivity), Internal = sqrt(diseases_curr$internal_connectivity), External = sqrt(diseases_curr$external_connectivity))
#     if (any(diseases_curr$Internal !=0)) {
#       j = j+1
#       diseases = bind_rows(diseases, diseases_curr)
# 
#     }
#     }
# 
#       mod_rand = lm(External ~ Internal-1, data = diseases)
#       rand_slope = coef(mod_rand)[[1]]
# 
#       coherence = (dis$Slope - rand_slope)/(median(msigs$Slope, na.rm = TRUE)-rand_slope)
#       coherence = data.frame(Disease = dis$Disease, Coherence = coherence)
# 
#       coherence_tab = bind_rows(coherence, coherence_tab)
#       saveRDS(coherence_tab, "./data/Coherence_Results/coherence_string_filt_CC.rds")
# }
# 

