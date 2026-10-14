rm(list = ls(all.names = TRUE)) # clears global environ.

# Load packages
library(dplyr) 
library(broom)
library(igraph)
library(tibble)
library(ggtext)
library(ggplot2)
library(tidyverse)
library(GSEABase)
library(clusterProfiler)

`%ni%` <- negate(`%in%`)

load("./Input/GSEA_data.RData")
load("./Input/Network_data.RData")

# Reformat GSEA and graph results
GESE_results <- data.frame(GSEA_results[-1], row.names = GSEA_results$ID)
community_results <- data.frame(Term = V(graph_subset)$name, Cluster = fine_grained_communities$membership)

# For each cluster, calculate penalized relationships between terms and ancestors.
Cluster_enrichment <- lapply(unique(community_results$Cluster), function(cluster){
  
  # Isolate GSEA enriched terms in current cluster
  terms <- community_results[which(community_results$Cluster == cluster), 'Term']
  terms <- terms[which(terms %in% GSEA_results$ID)]

  # It is possible that no GSEA enriched terms are found within a cluster
  # If none are found, return empty.
  if(length(terms) == 0){return()}
  
  all_nodes <- list()
  
  for(term in terms){
    
    ancestors <- unlist(go_obo$ancestors[term])
    all_nodes[[length(all_nodes) + 1]] <- ancestors[ancestors != term]
    
  }
  
  # Define ancestors for each term in the current cluster
  all_nodes <- unique(unlist(all_nodes))
  
  # Create an empty matrix of predefined size with GSEA terms in 
  # the current cluster as rows and ancestors as columns
  dist_mat <- matrix(NA, length(terms), length(all_nodes))
  rownames(dist_mat) <- terms
  colnames(dist_mat) <- all_nodes
  
  for(i in 1:nrow(dist_mat)){
    for(j in 1:ncol(dist_mat)){
      
      # Calculate the distance between each GSEA term and each ancestor term
      distance <- all_simple_paths(graph_subset,
                                   from = rownames(dist_mat)[i],
                                   to = colnames(dist_mat)[j],
                                   mode = "in")
      
      # If no direct upstream paths are found, distance is 0.
      # If an upstream path is found, subtract 1 as the search term is included in the list.
      distance_value <- ifelse(length(distance) > 0,
                               min(lengths(distance)) - 1,
                               0)
      
      # Isolate the IC score and GSEA term NES (Normalized Enrichment Score) for the current term
      IC_score <- go_ids_subset[which(go_ids_subset$go_id == colnames(dist_mat)[j]), 'ic_scores'][1]
      Term_NES <- abs(GESE_results[rownames(dist_mat)[i], 'NES'])
      
      # Calculate the weighted distance score between the search term and the current ancestor.
      # IC score (high means more specific) and term NES (normalized enrichment) are multiplied together, 
      # then scaled by Gaussian kernel decay of the distance between terms - 1. 
      # Distance - 1 was selected to ensure a distance of 1 (direct child-parent relation) is given a scaling weight of 1.
      dist_mat[i, j] <- ifelse(distance_value != 0,
                               IC_score * Term_NES * exp(-((distance_value - 1)**2)/(2 * (1**2))), # Gaussian kernel decay
                               Inf)
      
    }
  }
  
  # If the value is infinite, no direct path was found. Set it to 0 to give no impact on weighting sum.
  dist_mat[is.infinite(dist_mat)] <- 0
  
  # Sum weights by column. Higher weight mean more/closer connections between terms.
  sorted_list <- data.frame(colSums(dist_mat)) %>%
    rownames_to_column() %>%
    setNames(c('go_id', 'Score'))
  
  # Merge with GO term scoring metrics, arrange by score, and isolate top 5 terms.
  sorted_list <- merge(sorted_list,
                       go_ids_subset[!duplicated(go_ids_subset$go_id),
                                 c("go_id", "n_children", "ic_scores", "min_depth", "max_depth", "name_1006")],
                       all.x = TRUE,
                       by = 'go_id') %>%
    arrange(desc(Score)) %>% 
    slice_head(n = 5)
  
  # Add cluster details
  sorted_list$Cluster <- cluster
  
  # Print notice that cluster has been completed
  print(paste('Completed cluster', cluster, 'of', length(unique(community_results$Cluster))))
  
  return(sorted_list)
  
})

enriched_terms <- do.call(bind_rows, Cluster_enrichment)

# An enriched term can be enriched in more than one cluster. In that case, keep one.
enriched_terms <- enriched_terms %>% distinct(name_1006, .keep_all = TRUE)

########################################

tab_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
                "#bcbd22", "#17becf")

names(tab_colors) <- seq(1, length(tab_colors), by = 1)

# Limit GO terms to enriched terms
subset_go_ids <- go_ids_subset[which(go_ids_subset$go_id %in% enriched_terms$go_id), ]
subset_gsea_go_ids <- subset_go_ids[, c('go_id', 'hgnc_symbol')] %>% setNames(c('term', 'gene'))
subset_gsea_go_names <- subset_go_ids[, c('go_id', 'name_1006')] %>% setNames(c('term', 'name'))

########################################

# Run fresh GSEA to test enriched terms
gsea_results <- GSEA(geneList = gene_list,
                     TERM2GENE = subset_gsea_go_ids,
                     TERM2NAME = subset_gsea_go_names,
                     minGSSize = 1,            
                     maxGSSize = 5000,
                     pvalueCutoff  = 1,
                     pAdjustMethod = "BH",
                     verbose = FALSE)

########################################

# Determine which genes in the gene set were also DEGs
gene_sets <- gsea_results@geneSets

sig_in <- numeric(nrow(gsea_results@result))

for(i in 1:nrow(gsea_results@result)){
  
  sig_in[[i]] <- length(gene_sets[gsea_results@result$ID[i]][unlist(gene_sets[gsea_results@result$ID[i]]) %in% Sig_degs])
  
}

gsea_results@result$sig_in <- sig_in

########################################

# Add GO metrics to results
gsea_results@result <- merge(gsea_results@result, 
                             go_ids_subset[!duplicated(go_ids_subset$go_id), 
                                       c("go_id", "n_children", "ic_scores", "min_depth", "max_depth")] %>% 
                               setNames(c("ID", "n_children", "ic_scores", "min_depth", "max_depth")), 
                             by = "ID")

gsea_results@result$tags <- as.numeric(str_extract(gsea_results@result$leading_edge, "(?<=tags=).*?(?=%)"))
gsea_results@result$list <- as.numeric(str_extract(gsea_results@result$leading_edge, "(?<=list=).*?(?=%)"))
gsea_results@result$signal <- as.numeric(str_extract(gsea_results@result$leading_edge, "(?<=signal=).*?(?=%)"))
gsea_results@result$max_signal <- as.numeric(100)

########################################

# Prepare results for plotting
result <- gsea_results@result
enriched_terms <- enriched_terms[, c('name_1006', 'Cluster')] %>% setNames(c('Description', 'Cluster'))
result <- merge(result, enriched_terms, by = 'Description')
result$Description <- toupper(result$Description)
result <- result %>% arrange(Cluster, desc(NES))

# Create colored labels based on enrichment cluster
result$colored_label <- paste0("<span style='color:", 
                               tab_colors[as.character(result$Cluster)], 
                               ";'>", 
                               result$Description, 
                               "</span>")

result$colored_label <- factor(result$colored_label, levels = result$colored_label)

########################################

# Plot enriched GSEA data
plot <- ggplot(result, aes(x = NES, y = colored_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + # Line at x = 0
  geom_point(aes(size = max_signal), shape = 21, color = "grey50", show.legend = TRUE) + # Circle for 100% signal
  geom_point(aes(size = signal, color = max_depth), show.legend = FALSE) + # Circle for calculated signal and term depth
  
  theme_bw() +
  
  labs(size = 'Signal Strength') +
  
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.y.right = element_markdown(family = "mono", face = "bold", color = "black", size = 12),
        axis.text.y.left = element_markdown(family = "mono", face = "bold", color = "black", size = 12),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "right") +
  
  scale_y_discrete(labels = function(x){
    
    # Pad text with white space so graphs look good together
    padded <- str_pad(x,
                      width = max(str_width(go_ids_subset$name)),
                      side = "left",
                      pad = "_")
    
    sub("(_{2,})(.*?)(_{2,})", "<span style='color:#FFFFFF;'>\\1\\2\\3</span>", padded)
    
  },
  
  position = "left") +
  
  scale_color_gradient(low = "blue",
                       high = "red",
                       name = "Term Depth (Max)",
                       guide = guide_colorbar(frame.colour = "black", ticks = TRUE)) +
  
  scale_size_continuous(limits = c(0, 100), range = c(1, 10))

######

# Add lines seperating enrichment clusters
intercept = 0

for(i in 1:(max(result$Cluster) - 1)){
  
  intercept = intercept + nrow(result[which(result$Cluster == i), ])
  
  plot <- plot +
    geom_hline(yintercept = intercept + 0.5, linetype = "dashed", color = "black", linewidth = 1)
  
}

######

if(!dir.exists(paste0('./Output/Enriched_Graph_plots/'))) {
  dir.create(paste0('./Output/Enriched_Graph_plots/'), recursive = TRUE)
}

ggsave(paste0('./Output/Enriched_Graph_plots/subset_gsea_dotplot.png'),
       plot = plot,
       width = 20,
       height = 10,
       dpi = 300)
