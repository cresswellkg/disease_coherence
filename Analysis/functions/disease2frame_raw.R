diseases2frame_raw = function(disease_dis, links, disease_dis_name) {
  require(igraph)
  source('functions/genes2interactions.R')

  print(disease_dis)
  #Get internal vs. external connectivity
  
  disease_nets = lapply(disease_dis, function(x) {
    gene_list = disease_dis
    genes2interactions(gene_list, links, directed = FALSE)
  })
  
  
  disease_length = 1:length(disease_nets)
  
  #Get diseases
  
  disease_frame = lapply(disease_length, function(y) {
    x = disease_nets[[y]]
    sub_frame = bind_rows(x$internal_connectivity, x$external_connectivity)
    sub_frame = cbind.data.frame(colnames(sub_frame), t(sub_frame[1,]), t(sub_frame[2,]))
    colnames(sub_frame) = c("Genes", "Internal", "External")
    sub_frame %>% mutate(Disease = disease_dis_name[y])
  })

  disease_frame = bind_rows(disease_frame) %>% filter(!is.na(Disease)) 
  return(as.matrix(disease_frame))
}

