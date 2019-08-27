# remotes::install_github("schochastics/graphlayouts")
require(graphlayouts)
require(igraph)
require(ggraph)

#Creating ring networks
#10 vertices
ring1 = make_ring(10)
#20 vertices
ring2 = make_ring(20)

#Plotting each relationship
set.seed(1)
#To change line color change the "colour" parameter in geom_edge_link
#To change edge color change value in scale_fill_manual
ring1_plot = ggraph(ring1, layout = "stress") +
  geom_edge_link(width=.5,colour="dodgerblue", alpha = 0.5)+
  geom_node_point(shape=21,size=3, aes(fill = "green")) +
  theme_graph() + 
  scale_fill_manual(values = "darkorange") + 
  guides(fill = FALSE)

ring2_plot = ggraph(ring2, layout = "stress") +
  geom_edge_link(width=.5,colour="dodgerblue", alpha = 0.5)+
  geom_node_point(shape=21,size=3, aes(fill = "green")) +
  theme_graph() + 
  scale_fill_manual(values = "darkorange") + 
  guides(fill = FALSE)

#Full graphs

full1_plot = ggraph(make_full_graph(10), layout = "stress") +
  geom_edge_link(width=.5,colour="dodgerblue", alpha = 0.5)+
  geom_node_point(shape=21,size=3, aes(fill = "green")) +
  theme_graph() + 
  scale_fill_manual(values = "darkorange") + 
  guides(fill = FALSE)

full2_plot = ggraph(make_full_graph(20), layout = "gem") +
  geom_edge_link(width=.5,colour="dodgerblue", alpha = 0.5)+
  geom_node_point(shape=21,size=3, aes(fill = "green")) +
  theme_graph() + 
  scale_fill_manual(values = "darkorange") + 
  guides(fill = FALSE)

#Saving plots
ggsave(ring1_plot, file = "ring1.png", width = 2.5, height = 2.5, units = "in")
ggsave(ring2_plot, file = "ring2.png", width = 2.5, height = 2.5, units = "in")
ggsave(full1_plot, file = "full1.png", width = 2.5, height = 2.5, units = "in")
ggsave(full2_plot, file = "full2.png", width = 2.5, height = 2.5, units = "in")

