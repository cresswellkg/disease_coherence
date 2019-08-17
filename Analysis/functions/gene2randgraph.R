gene2randgraph = function(gene_list, actions, directed) {
  act_sub = actions %>% 
    filter(hgnc_symbol_a %in% gene_list | hgnc_symbol_b %in% gene_list) %>% mutate(hgnc_acting = ifelse(a_is_acting == "t",hgnc_symbol_a, hgnc_symbol_b),
                                                                                   hgnc_receiving = ifelse(a_is_acting == "t",
                                                                                                           hgnc_symbol_b, hgnc_symbol_a)) %>% dplyr::select(hgnc_acting, hgnc_receiving) %>% filter(complete.cases(.))
  gene_graph = graph_from_data_frame(act_sub, directed = directed)
  rand_graph = sample_gnm(length(V(gene_graph)), length(E(gene_graph)), directed = directed)
  
  V(rand_graph)$name = V(gene_graph)$name
  
  sub_graph = induced_subgraph(rand_graph,vids = which(V(gene_graph)$name %in% gene_list))
  internal_con = degree(sub_graph, mode = "in", loops = FALSE)
  external_con = degree(rand_graph, mode = "out", v = V(rand_graph)[V(rand_graph)$name%in% gene_list], loops = FALSE)
  return(list(internal_connectivity = internal_con,external_connectivity = external_con))
  
  }