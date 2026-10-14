rm(list = ls(all.names = TRUE)) # clears global environ.

# Load packages
library(dplyr) 
library(broom)
library(tibble)
library(igraph)
library(tidyverse)
library(org.Hs.eg.db)
library(ontologyIndex)
library(ontologySimilarity)

# Create not in function by negating in
`%ni%` <- negate(`%in%`)

# Check if GO to gene symbol mapping has been complete and is saved in ./Input
# If not, map gene symbols to related GO terms, specifying GOALL, which maps all
# genes related to a given GO term including those linked to upstream ancestors.
if(file.exists("./Input/GOALL.RData") == FALSE){
  
  # Get all Entrez Gene IDs available in the database
  gene_universe <- data$HUGO_IDs
  
  end <- 1
  query <- list()
  for(i in 1:ceiling(length(gene_universe)/500)){
    
    start <- end
    end <- 500*i
    
    # Retrieve GO-to-gene mappings
    go_map <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = gene_universe[start:end],
                                    columns = c("SYMBOL", "GOALL"),
                                    keytype = "SYMBOL")
    
    go_map <- distinct(go_map, paste(SYMBOL, GOALL), .keep_all = TRUE)
    query[[i]] <- go_map[, c("SYMBOL", "GOALL")]
    
  }
  
  go_ids <- do.call(bind_rows, query)
  
  save(go_ids,
       gene_universe,
       file = './Input/GOALL.RData')
  
} else{load('./Input/GOALL.RData')}

################################################################################
# Load differential expression data. The example in this repository is a modified version of data available through:

# "Leveraging a comprehensive unbiased RNAseq database to characterize human monocyte-derived 
# macrophage gene expression profiles within commonly employed in vitro polarization methods" 
# published in Scientific Reports at: https://doi.org/10.1038/s41598-024-78000-6 

# Project code is available at: https://github.com/trsmyth/Comprehensive-Unbiased-RNAseq-Database-Uncovers-New-Human-Monocyte-Derived-Macrophage-Subtypes

# Specifically, 3. Differential Expression/Figure 1-2 and Supplemental Figure S3-4 - Differential Expression by Polarization Method.Rmd
# Differential expression was conducted as described, but all output columns were retained and gene symbols were merged with DE results.
# This provides access to the moderated t-statistic which is used for GSEA analysis.
################################################################################

load("./Input/DEG_result_example.RData")
data <- M_LPS_IFN_vs_M_LPS %>% distinct(HUGO_IDs, .keep_all = TRUE)

########################################

# GO term ancestor-child relationships are mapped based. 
# Requires .obo file from https://geneontology.org/docs/download-ontology/ 
go_obo <- get_ontology(file = './Input/go-basic.obo', extract_tags = "everything")

# Information content value calculated using term frequencies
frequencies <- get_term_frequencies(go_obo, go_obo[["id"]])
ic_scores <- get_term_info_content(go_obo, names(frequencies))

########################################

# Create edge list linking parents
edges_list <- lapply(names(go_obo$parents), function(term){
  
  # ID parents for current term
  parents <- go_obo$parents[[term]]
  
  # If there is at least one parent, track from current term to every parent
  if(length(parents) > 0){data.frame(from = term, to = parents, stringsAsFactors = FALSE)} 
  else{NULL}
  
})

edges_df <- do.call(rbind, edges_list) # Combine resulting term to parent relationships
graph <- graph_from_data_frame(d = edges_df, directed = TRUE) # Create a graph linking these relationships up the tree

# Find root nodes (BP, CC, MF)
all_nodes <- V(graph)$name
out_degrees <- igraph::degree(graph, mode = "out")
roots <- names(out_degrees[out_degrees == 0])[1:3]

########################################

# Define gene universe using all genes tested for DE
gene_universe <- rownames(data) 

# Filter GO terms requiring at least 5 genes per term
go_ids <- go_ids[, c('SYMBOL', 'GOALL')] %>% setNames(c('hgnc_symbol', 'go_id'))
go_ids <- go_ids[which(go_ids$go_id != '' & go_ids$hgnc_symbol != ''), ]
go_ids <- go_ids %>% group_by(go_id) %>% filter(n()>= 5)

# Determine key graph information, including term depth, child count, and IC score.
graph_data <- lapply(unique(go_ids$go_id), function(single_go_id){
  
  # Determine the number of children terms the current term has
  n_children <- length(go_obo[['children']][[single_go_id]])
  
  # Isolate the IC (information content) score
  ic_scores <- ic_scores[single_go_id]
  
  # Determine if term is obsolete
  obsolete <- as.character(go_obo[["obsolete"]][single_go_id])
  
  #########
  
  path_depth <- tryCatch({
    
    # Find reachable roots (BP, CC, MF)
    reachable_roots <- intersect(roots, subcomponent(graph, single_go_id, mode = "out")$name)
    
    # Calculate all paths to each reachable root
    all_paths <- lapply(reachable_roots, function(root){
      
      # Maps non-repeating route from node to roots
      all_simple_paths(graph, from = single_go_id, to = root, mode = "out")
      
    })
    
    all_paths <- unlist(all_paths, recursive = FALSE)
    
    # Calculate each path length - 1
    path_lengths <- sapply(all_paths, length) - 1 
    min_path_depth <- min(path_lengths)
    max_path_depth <- max(path_lengths)
    
    list(go_obo[["name"]][reachable_roots], min_path_depth, max_path_depth)
    
  }, 
  error = function(e){
    
    reachable_roots <- intersect(roots, subcomponent(graph, single_go_id, mode = "out")$name)
    min_path_depth <- NA
    max_path_depth <- NA
    
    list(go_obo[["name"]][reachable_roots], min_path_depth, max_path_depth)
    
  })
  
  #########
  
  results <- data.frame('go_id' = single_go_id, 
                        'ontology' = path_depth[[1]],
                        'n_children' = n_children, 
                        'ic_scores' = ic_scores, 
                        'obsolete' = obsolete, 
                        'min_depth' = path_depth[[2]], 
                        'max_depth' = path_depth[[3]])
  
  return(results)
  
})

# Combine GO terms with key graph information
graph_data <- do.call(bind_rows, graph_data)
go_ids <- merge(go_ids, graph_data, by = 'go_id')
go_ids$name_1006 <- go_obo[["name"]][go_ids$go_id]

########################################

# Remove obsolete terms and terms without a description
go_ids_subset <- go_ids[which(go_ids$obsolete == 'FALSE'), ]
go_ids_subset <- go_ids_subset[which(go_ids_subset$name_1006 != ''), ]

go_ids_subset <- go_ids_subset[which(go_ids_subset$ontology == 'biological_process'), ] # Isolate BP terms
gsea_go_ids <- go_ids_subset[, c('go_id', 'hgnc_symbol')] %>% setNames(c('term', 'gene')) # Get ID to gene relations
gsea_go_names <- go_ids_subset[, c('go_id', 'name_1006')] %>% setNames(c('term', 'name')) # Get ID to term name relations

########################################

save(file = './Input/GO_setup.RData', 
     data,
     go_obo,
     go_ids_subset, 
     gsea_go_ids, 
     gsea_go_names)
