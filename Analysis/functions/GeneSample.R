GeneSample = function(num_genes, actions, directed = TRUE) {
  
  require(igraph)
  gene_list = sample(unique(actions$hgnc_symbol_a), num_genes)
  
  if (directed) {
  act_sub = actions %>% 
    filter(hgnc_symbol_a %in% gene_list | hgnc_symbol_b %in% gene_list) %>% mutate(hgnc_acting = ifelse(a_is_acting == "t",hgnc_symbol_a, hgnc_symbol_b),
                                                                                   hgnc_receiving = ifelse(a_is_acting == "t",
                                                                                                          hgnc_symbol_b, hgnc_symbol_a)) %>% dplyr::select(hgnc_acting, hgnc_receiving) %>% filter(complete.cases(.))
  gene_graph = graph_from_data_frame(act_sub, directed = directed)
  sub_graph = induced_subgraph(gene_graph,vids = which(V(gene_graph)$name %in% gene_list))
  internal_con = degree(sub_graph, mode = "in", loops = FALSE)
  external_con = degree(gene_graph, mode = "out", v = V(gene_graph)[V(gene_graph)$name%in% gene_list], loops = FALSE)
  }
  
  else {
    link_sub = actions %>% 
      filter(hgnc_symbol_a %in% gene_list | hgnc_symbol_b %in% gene_list) %>% 
      mutate(hgnc_symbol_a = ifelse(hgnc_symbol_a == "", NA, hgnc_symbol_a),
             hgnc_symbol_b = ifelse(hgnc_symbol_b == "", NA, hgnc_symbol_b)) %>%
      filter(complete.cases(.))  %>% dplyr::select(hgnc_symbol_a, hgnc_symbol_b) %>%
      mutate(key = paste0(pmin(hgnc_symbol_a, hgnc_symbol_b), pmax(hgnc_symbol_a, hgnc_symbol_b), sep = "")) %>% 
      group_by(key) %>% dplyr::slice(1) %>% ungroup() %>% dplyr::select(-key)
    
    gene_graph = graph_from_data_frame(link_sub, directed = directed)
    sub_graph = induced_subgraph(gene_graph,vids = which(V(gene_graph)$name %in% gene_list))
    internal_con = degree(sub_graph, loops = FALSE)
    external_con = degree(gene_graph, v = V(gene_graph)[V(gene_graph)$name%in% gene_list], loops = FALSE)
    
  }
  
  return(list(graph = gene_graph, genes = gene_list, internal_connectivity = internal_con,external_connectivity = external_con))
  
}
