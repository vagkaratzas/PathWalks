# PathWalks
PathWalks is a random walk based algorithm where a walker crosses a pathway-to-pathway network under the guidance of a disease-related map. The latter is a gene network that we construct by integrating multi-source information regarding a specific disease. The most frequent trajectories highlight communities of pathways that are expected to be strongly related to the disease under study.

Required library: igraph
Arguments: gene_network_file, pathway_network_file, converging_factor, gene_restart_timer, check_similarity_timer, number_of_last_variances_to_check, converging_variance
Network files must be in edgelist form.
Example command line call: Rscript pathwalks.R "geneEdgelistIPF.tsv" "hsaPathwayEdgelistIPF.tsv" 0.95 50 100 10 0.003