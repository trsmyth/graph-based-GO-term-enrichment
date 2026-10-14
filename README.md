# graph-based-GO-term-enrichment
Graph-based approach for the identification and collapse of related GO terms following enrichment analysis

Pathway enrichment techniques are commonly employed following differential expression (DE) analysis to identify biological themes within groups of interest. However, common techniques like gene ontology (GO) enrichment or gene set enrichment analysis (GSEA) rely on a priori pathway selection or extensive post-enrichment filtering to remove irrelevant pathways. For example, airway epithelial cells frequently demonstrate enrichment for pathways related to sperm motility due to shared expression of ciliary genes while blood-derived macrophages demonstrate enrichment of T-cell pathways, likely due to the extensive volume of research linked to T-cells over macrophages.

This repository contains code related to a graph-based enrichment approach utilizing the GO Directed Acyclic Graph (DAG). GO parent-child relationships can be defined, nodes (terms) can be clustered, and enriched themes within specific clusters can be determined to identify latent themes within enrichment. The approach described here utilizes differential expression results to conduct GSEA of GO pathways which allows for calculation of a weighted latent enrichment, however, this technique would also work with pure GO enrichment or other graph-based enrichment technique.

Below is an example enrichment workflow. Data used for enrichment is from "Leveraging a comprehensive unbiased RNAseq database to characterize human monocyte-derived macrophage gene expression profiles within commonly employed in vitro polarization methods" published in Scientific Reports at: https://doi.org/10.1038/s41598-024-78000-6, and the specific input data is available in the 'Input' folder.

Specifically, this data are DE results from macrophages polarized with interferon gamma and LPS versus LPS alone. Initial GSEA using the moderated t-statistic for gene weighting demonstrates a large number of enriched pathways (limited to 50 pathways). However, many highly related terms and seemingly irrelevant terms, such as T-cell linked pathways, are enriched.
<br><br>
<p align="center">
  <img width="1000" alt="image" src="https://github.com/user-attachments/assets/9e5d2dcc-467c-42c5-95b2-3c144ecd71fa" />
</p>

A new graph is formed which links enriched terms to all known ancestors, producing clusters of enrichment which can be assessed for latent themes.  
  <br>
<p align="center">
  <img width="500" alt="image" src="https://github.com/user-attachments/assets/9a297476-c325-4d67-8b75-e1f0bdfe4fc9" />
</p>

Terms within each cluster are isolated and the step distance between significant GSEA terms and all possible ancestors are assessed. Ancestor term enrichment is assessed using a weighted score with the following formula:  

<div align="center">
  \$`\text{IC score} \times \text{NES} \times \exp\left(-\frac{(\text{distance} - 1)^2}{2\sigma^2}\right)`\$
</div>
<br>
- IC score = GSEA term information content score (higher is a more specific term)
<br>
- NES = normalized enrichment score
<br>
- Gaussian kernel decay is utilized where distance is the step distance between the significant GSEA term and the current ancestor. Sigma is set to 1 and distance is subtracted by 1 to ensure direct parent-child relationships have a weight of 1.
<br><br>
Weighted scores are assessed for every significant GSEA term within the cluster with each potential ancestor within the cluster (where no connections are set to 0) and enriched ancestor terms are selected based on the top n summed scores. The selected ancestor terms can then be assessed by GSEA as well and graphed based on cluster to identify latent themes within the data.
<br><br>
<p align="center">
  <img width="1000" alt="image" src="https://github.com/user-attachments/assets/e809bc13-456d-4b10-b4e5-d7fc592efd14" />
</p>
