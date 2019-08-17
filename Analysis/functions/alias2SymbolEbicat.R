#data(ebicat37)
#Function for converting ebicat37 mapped genes data to symbols designed specifically for the ebicatdata

alias2SymbolEbicat = function(ebicat37) {
  
  ebi_dat = as.data.frame(mcols(ebicat37))
  
  #Fix spaces between - and ,
  
  ebi_dat = unembed(ebi_dat, "MAPPED_GENE", sep = ",")
  
  ebi_dat$MAPPED_GENE = gsub(" -", "-", ebi_dat$MAPPED_GENE)
  ebi_dat$MAPPED_GENE = gsub("- ", "-", ebi_dat$MAPPED_GENE)
  
  #Summarise the results
  
  #Keeping track of changes (sapply necessary to keep placeholders for unmapped genes)
  
  recomb_genes = sapply(ebi_dat$MAPPED_GENE, function(x) alias2Symbol(x, species = "Hs"))
  
  recomb_genes = unlist(lapply(recomb_genes, function(x) paste0(x,collapse=' ')))
  
  change = (ebi_dat$MAPPED_GENE == recomb_genes)
  
  sum_vec = ifelse(change == TRUE, "Same", ifelse(change == FALSE & recomb_genes == "", "No Match", "Converted"))
  
  ebi_dat$MAPPED_GENE = recomb_genes
  
  return(list(ebi_dat = ebi_dat, change = change))
}