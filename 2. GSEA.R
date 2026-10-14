rm(list = ls(all.names = TRUE)) # clears global environ.

# Load packages
library(dplyr) 
library(broom)
library(tibble)
library(ggtext)
library(ggplot2)
library(GSEABase)
library(tidyverse)
library(clusterProfiler)

`%ni%` <- negate(`%in%`)

########################################

load('./Input/GO_setup.RData')

n_pathways = 50

# Isolate significant DEGs
Sig_degs <- data[which(data$adj.P.Val < 0.05), 'HUGO_IDs']

# Arrange DE results by descending modified t-statistic
gene_list <- data %>% 
  arrange(desc(t)) %>% 
  pull(t, name = HUGO_IDs)

#########

# Run GSEA. An alternative for this workflow would be genering GO enrichment from a list of genes.
gsea_results <- GSEA(geneList = gene_list,
                     TERM2GENE = gsea_go_ids,
                     TERM2NAME = gsea_go_names,
                     minGSSize = 5,            
                     maxGSSize = 100,
                     pvalueCutoff  = 0.05,
                     pAdjustMethod = "BH",
                     verbose = FALSE)

#########

# Determine the number of DE genes are in each defined gene set, i.e. GO term.
gene_sets <- gsea_results@geneSets

sig_in <- numeric(nrow(gsea_results@result))

for(i in 1:nrow(gsea_results@result)){
  
  sig_in[[i]] <- length(gene_sets[gsea_results@result$ID[i]][unlist(gene_sets[gsea_results@result$ID[i]]) %in% Sig_degs])
  
}

gsea_results@result$sig_in <- sig_in

#########

# Merge GSEA results and key graph information from GO data.
gsea_results@result <- merge(gsea_results@result, 
                             go_ids_subset[!duplicated(go_ids_subset$go_id), 
                                           c("go_id", "n_children", "ic_scores", "min_depth", "max_depth")] %>% 
                               setNames(c("ID", "n_children", "ic_scores", "min_depth", "max_depth")), 
                             by = "ID")

gsea_results@result$tags <- as.numeric(str_extract(gsea_results@result$leading_edge, "(?<=tags=).*?(?=%)"))
gsea_results@result$list <- as.numeric(str_extract(gsea_results@result$leading_edge, "(?<=list=).*?(?=%)"))
gsea_results@result$signal <- as.numeric(str_extract(gsea_results@result$leading_edge, "(?<=signal=).*?(?=%)"))
gsea_results@result$max_signal <- as.numeric(100)

#########

# Arrange by descending max depth (how many branches down, i.e. how many parents) and isolate top n_pathways
result <- gsea_results@result %>%
  arrange(desc(max_depth)) %>%
  slice_head(n = n_pathways) %>%
  arrange(desc(NES))

result$Description <- toupper(result$Description)

#########

# Plot GSEA results
plot <- ggplot(result, aes(x = NES, y = Description)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + # Set vline to separate neg and pos pathways
  geom_point(aes(size = max_signal), shape = 21, color = "grey50", show.legend = TRUE) + # Create outer point for max signal (100%)
  geom_point(aes(size = signal, color = max_depth), show.legend = FALSE) + # Create inner point for calculated signal colored by pathway depth
  
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
    
    # Pad pathway names to ensure plots are uniform in size if multiple enrichment are run.
    # This is solely for presentation purposes to make related graphs look good in the same panel.
    padded <- str_pad(x,
                      width = max(str_width(go_ids_subset$name_1006)),
                      side = "left",
                      pad = "_")
    
    sub("(_{2,})(.*?)(_{2,})", "<span style='color:#FFFFFF;'>\\1\\2\\3</span>", padded)
    
  },
  
  position = "left") +
  
  # Create color gradient scale for term depth data
  scale_color_gradient(low = "blue",
                       high = "red",
                       name = "Term Depth (Max)",
                       guide = guide_colorbar(frame.colour = "black", ticks = TRUE)) +
  
  # Create continuous scale from 100% to 0% for data points
  scale_size_continuous(limits = c(0, 100), range = c(1, 10))

#########

# Save results, creating directory if not present
if(!dir.exists(paste0('./Output/GSEA_plots/'))){
  dir.create(paste0('./Output/GSEA_plots/'), recursive = TRUE)
  }

ggsave(paste0('./Output/GSEA_plots/', 'gsea_dotplot.png'),
       plot = plot,
       width = 20,
       height = 10,
       dpi = 300)

#########

# For exporting results
GSEA_results <- result

save(file = './Input/GSEA_data.RData', 
     GSEA_results, 
     go_ids_subset, 
     go_obo, 
     Sig_degs, 
     gene_list)
