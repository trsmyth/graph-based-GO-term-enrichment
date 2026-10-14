rm(list = ls(all.names = TRUE)) # clears global environ.

# Load packages
library(dplyr) 
library(broom)
library(igraph)
library(tibble)
library(ggplot2)
library(tidyverse)

load('./Input/GSEA_data.RData')

# Isolate GO pathway ID, name, and linked gene symbols
go_ids_subset <- go_ids_subset[, c('go_id', 'name_1006', 'hgnc_symbol')] %>% setNames(c('term', 'name', 'gene'))

plot = TRUE

# Define all possible nodes including enriched terms and all their ancestors
all_nodes <- unique(unlist(go_obo$ancestors[GSEA_results$ID]))

# Create empty graph and add each node as a vertex
graph_subset <- make_empty_graph()
graph_subset <- add_vertices(graph_subset, nv = length(all_nodes), name = all_nodes)

# For every node, draw the edges between node and its ancestors
for(ID in all_nodes){
  
  term_parents <- unlist(go_obo$parents[ID])
  
  if(length(term_parents) > 0){
    for(parent in term_parents){
      if(parent %in% all_nodes){
        graph_subset <- add_edges(graph_subset, c(parent, ID))
      }
    }
  }
}

# Define labels of graph vertices. If term is enriched from GSEA, change to red. Set root (biological process) to green. 
V(graph_subset)$label.color <- ifelse(V(graph_subset)$name %in% GSEA_results$ID, "tomato", "grey50")
V(graph_subset)['GO:0008150']$label.color <- "forestgreen"

# Create the graph layout and conduct louvain clustering
layout = layout_with_fr(graph_subset)
fine_grained_communities <- cluster_louvain(as_undirected(graph_subset), resolution = 0.25)

# Plot if desired
if(plot == TRUE){
  
  if(!dir.exists(paste0('./Output/Graph_plots/'))) {
    dir.create(paste0('./Output/Graph_plots/'), recursive = TRUE)
  }
  
  png(file = paste0('./Output/Graph_plots/louvain_plot.png'), height = 5000, width = 7500, res = 300)
  
  plot(fine_grained_communities,
       graph_subset,
       main = "Louvain Community Detection",
       vertex.label.cex = 0.5,
       edge.label.cex = 0.5,
       layout = layout,
       vertex.label.angle = -pi*10,
       vertex.size = 0.25,
       edge.arrow.size = 0,
       margin = -0.05,
       vertex.label.color = graph_subset$label.color)
  
  dev.off()
  
  png(file = paste0('./Output/Graph_plots/graph_plot.png'), height = 5000, width = 7500, res = 300)
  
  plot(graph_subset,
       vertex.label.cex = 0.5,
       edge.label.cex = 0.5,
       layout = layout,
       vertex.label.angle = -pi*10,
       vertex.size = 0.25,
       edge.arrow.size = 0,
       margin = -0.1,
       vertex.label.color = graph_subset$label.color)
  
  dev.off()
  
}

save(file = './Input/Network_data.RData', graph_subset, fine_grained_communities)
