#Function for converting ebicat37 mapped genes data to symbols designed specifically for the ebicatdata
#You must first run data(ebicat37) from the gwascat package


alias2SymbolEbicat = function(ebicat37) {
  
  #Convert ebicat37 to a data frame
  ebi_dat = as.data.frame(mcols(ebicat37))
  
  #Fix spaces between - and ,
  
  ebi_dat = unembed(ebi_dat, "MAPPED_GENE", sep = ",")
  
  ebi_dat$MAPPED_GENE = gsub(" -", "-", ebi_dat$MAPPED_GENE)
  ebi_dat$MAPPED_GENE = gsub("- ", "-", ebi_dat$MAPPED_GENE)
  
  #Summarise the results
  
  #Keeping track of changes (sapply necessary to keep placeholders for unmapped genes)
  
  recomb_genes = sapply(ebi_dat$MAPPED_GENE, function(x) alias2Symbol(x, species = "Hs"))
  
  #Find where there are multiple genes mapped to a symbol and combine them with a space delimeter
  recomb_genes = unlist(lapply(recomb_genes, function(x) paste0(x,collapse=' ')))
  
  #Identify location of changed genes
  change = (ebi_dat$MAPPED_GENE == recomb_genes)
  
  #Create a vector containing all changed values
  sum_vec = ifelse(change == TRUE, "Same", ifelse(change == FALSE & recomb_genes == "", "No Match", "Converted"))
  
  #Replace old genes with the new mapped genes
  ebi_dat$MAPPED_GENE = recomb_genes
  
  return(list(ebi_dat = ebi_dat, change = change))
}