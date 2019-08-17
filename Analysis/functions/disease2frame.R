diseases2mod = function(disease_dis, links) {
  
  source('functions/genes2graph.R')

  
  #Get internal vs. external connectivity
  
  disease_nets = lapply(disease_dis, function(x) {
    gene_list = read.table(x)$V1
    genes2graph(gene_list, links, directed = FALSE)
  })
  
  
  disease_length = 1:length(disease_nets)
  
  #Get diseases
  
  disease_dis_name = gsub('.*/', '', disease_dis)
  
  disease_dis_name = gsub('.txt', '', disease_dis_name)
  
  
  
  disease_frame = lapply(disease_length, function(y) {
    x = disease_nets[[y]]
    sub_frame = bind_rows(x$internal_connectivity, x$external_connectivity)
    sub_frame = cbind.data.frame(colnames(sub_frame), t(sub_frame[1,]), t(sub_frame[2,]))
    colnames(sub_frame) = c("Genes", "Internal", "External")
    sub_frame %>% mutate(Disease = disease_dis_name[y])
  })
  
  disease_frame = bind_rows(disease_frame)
  return(disease_frame)
}

