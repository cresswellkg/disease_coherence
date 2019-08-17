diseases2mod = function(disease_dis, links) {
  require(igraph)
  source('functions/genes2graph.R')

  
  #Get internal vs. external connectivity
  
  disease_nets = lapply(disease_dis, function(x) {
    gene_list = read.table(x)$V1
    genes2graph(gene_list, links, directed = FALSE)
  })
  
  clusters = cluster_louvain(disease_nets[[1]])
  
  modularity(disease_nets[[1]], clusters$membership)
  
 
}

