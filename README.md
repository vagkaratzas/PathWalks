# PathWalks
PathWalks is a random walk based algorithm, developed in R, where a walker crosses a pathway-to-pathway network under the guidance of a disease-related map. The latter is a gene network that we construct by integrating multi-source information regarding a specific disease. The most frequent trajectories highlight communities of pathways that are expected to be strongly related to the disease under study.

- Required library: igraph
- Arguments: gene_network_file, pathway_network_file, converging_factor, gene_restart_timer, check_similarity_timer, number_of_last_variances_to_check, converging_variance
- Network files must be in edgelist form.
- Gene network must contain gene symbols.
- Pathway network must contain KEGG hsa pathway ids (example: hsaPathwayEdgelist_commonGenes.tsv, showing functional connections among pathways as found in KEGG and each edge weight represents the number of common genes between two pathways).
- Results: 1. Pathway ranks 2. Re-weighted pathways' network edgelist 3. Pathway clusters
- Example command line call: Rscript pathwalks.R "geneEdgelistIPF.tsv" "hsaPathwayEdgelistIPF.tsv" 0.95 50 100 10 0.003

The three translation files must be in the same directory with the pathwalks.R script.

Use case example input data included.
- Idiopathic Pulmonary Fibrosis (geneEdgelistIPF.tsv, hsaPathwayEdgelistIPF.tsv)
- Alzheimer's Disease (geneEdgelistAD.tsv, hsaPathwayEdgelistAD.tsv)

For more information read our article in https://www.biorxiv.org/content/10.1101/2020.01.27.921270v1.full.pdf.
