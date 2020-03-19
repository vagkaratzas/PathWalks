#!/usr/bin/env Rscript
# Random Walks on a Pathways' Network guided by Gene Network Information from  GeneMANIA results
# Required library: igraph
# Arguments: gene_network_file, pathway_network_file, converging_factor, gene_restart_timer, check_similarity_timer, number_of_last_variances_to_check, converging_variance
# Network files must be in edgelist form.
# Results: 1. Pathway ranks 2. Re-weighted pathways' network edgelist 3. Pathway clusters
# Example command line call: Rscript pathwalks.R "geneEdgelistIPF.tsv" "hsaPathwayEdgelistIPF.tsv" 0.95 50 100 10 0.003

translateInputGenes <- function(filename) {
  geneSymbolEdgelist <- as.matrix(read.delim(filename, header = FALSE))
  hsaGeneEdgelist <- matrix("", nrow = 0, ncol = 3) # initializing output translated hsa genes' matrix
  for (i in 1:nrow(geneSymbolEdgelist)){
    # must find both genes to translate, else removing edge
    gene1 <- geneSymbolEdgelist[i, 1][[1]]
    gene2 <- geneSymbolEdgelist[i, 2][[1]]
    hsaGene1 <- gene_translation_matrix[which(gene_translation_matrix[, 2] %in% gene1), 1][1] # in the case of a gene translating to more than one ids, we keep the first id
    hsaGene2 <- gene_translation_matrix[which(gene_translation_matrix[, 2] %in% gene2), 1][1]
    if (!is.na(hsaGene2) && !is.na(hsaGene1)) hsaGeneEdgelist <- rbind(hsaGeneEdgelist, c(as.character(hsaGene1[1]), as.character(hsaGene2[1]), geneSymbolEdgelist[i, 3]))
  }
  return(hsaGeneEdgelist)
}

# creating gene-map from the translated gene edgelist
createGeneNetwork <- function(edgelist) {
  graph <- graph_from_edgelist(edgelist[, 1:2], directed = FALSE)
  E(graph)$weight <- as.double(edgelist[, 3])
  # remove loops and multiple edges, simplify sum aggregates same edges
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "sum"))
  return(graph)
}

# creating pathways' network from input hsa KEGG pathways' edgelist
createPathwayNetwork <- function(edgelist) {
  pathGr <- graph_from_edgelist(edgelist[, 1:2], directed = FALSE)
  E(pathGr)$weight <- as.double(edgelist[, 3])
  # remove loops and multiple edges, simplify sum aggregates same edges
  pathGr <- simplify(pathGr, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "sum"))
  E(pathGr)$weight <- 1/E(pathGr)$weight # inverting to correctly calculate shortest paths
  return(pathGr)
}

# removing unconnected subgraphs
removeDetachedNodes <- function(graph){
  # i from 2, because we only keep the main graph
  for (i in 2:components(graph)$no){
    to_remove <- names(components(graph)$membership[components(graph)$membership==i])
    graph <- delete_vertices(graph, to_remove)
  }
  return(graph)
}

# levy flight with cauchy distribution random walk with restarts
geneStepMultipleRWR <- function(gG, gI, timer) {
  # cat(sprintf("GENE STEP RWR\n"))
  if (timer %% gene_restart_timer == 0){ # mod function for restart based on timer
    # cat(sprintf("Random Gene Restart\n"))
    gI <- sample(1:length(V(gG)), 1) 
  }
  x <- ceiling(abs(rcauchy(1))) # decide on number of map steps based on a cauchy distribution (levy flight)
  if (x > 10) x <- 10 # 10 steps max, to keep execution time low
  # cat(sprintf("Cauchy x = %d\n", x))
  gene_matrix <- V(gG)[gI]$name # first element of gene path
  # cat(sprintf("First element of gene path: %s\n", gene_matrix))
  for (i in 1:x){
    # cat(sprintf("Cauchy step: %d / %d\n", i, x))
    neighbors <- as.matrix(neighbors(gG, gI))
    sum <- sum(E(gG)[gI %--% as.numeric(V(gG)[rownames(neighbors)])]$weight) # sum of neighbors' weights
    # cat(sprintf("Sum of neighbors' weights: %f\n", sum))
    rand <- runif(1, min = 0, max = sum) # using random uniform distribution due to decimal values
    # cat(sprintf("Uniform rand: %f\n", rand))
    # monte carlo sampling
    current_sum <- 0
    for (j in 1:length(neighbors)){
      current_sum <- current_sum + E(gG)[gI %--% as.numeric(V(gG)[rownames(neighbors)[j]])]$weight
      # cat(sprintf("current_sum: %f\n", current_sum))
      if (rand <= current_sum){
        # cat(sprintf("current_sum: %f\n", current_sum))
        gene_matrix <- c(gene_matrix, rownames(neighbors)[j])
        # cat(sprintf("Adding element to gene path: %s\n", rownames(neighbors)[j]))
        gI <- as.numeric(V(gG)[rownames(neighbors)[j]])
        # cat(sprintf("new gI: %d\n", gI))
        break
      }
    }
  }
  return(list(gene = unique(gene_matrix), gene_it = gI))
}

# stepping to next pathway based on a matrix of genes and incrementing output scores
pathwayStepMultiple <- function(pG, g, pI) {
  # cat(sprintf("PATHWAY STEP with %s Gene(s)", g))
  nP <- length(V(pG))
  current_pI <- pI # save current pathway id
  countHitPathways <- matrix(0, nrow = nP) # initialize array to increment potential pathways based on gene hits
  rownames(countHitPathways) <- V(pG)$name # assigning pathway names
  # search gene matrix entries in gene_pathway_links_matrix to fetch pathway candidates
  for (i in 1:length(g)){
    candidate_pathways <- gene_pathway_links_matrix[which(gene_pathway_links_matrix[, 2] %in% g[i]), 1]
    if (length(candidate_pathways) > 0){
      candidate_pathways <- candidate_pathways[which(candidate_pathways %in% rownames(countHitPathways))] # removing hit pathways that are not in network
      countHitPathways[candidate_pathways, ] <- countHitPathways[candidate_pathways, ] + 1
    }
  }
  # normalization
  countHitPathways <- cbind(rownames(countHitPathways), countHitPathways[, 1])
  countHitPathways <- countHitPathways[order(countHitPathways[, 1]),]
  rownames(countHitPathways) <- c()
  countHitPathways <- cbind(countHitPathways, as.numeric(as.character(num_genes_per_pathway[, 2]))) # num_genes_per_pathway pre-sorted to match respective pathways
  countHitPathways <- cbind(countHitPathways, as.numeric(as.character(countHitPathways[, 2])) / as.numeric(as.character(countHitPathways[, 3]))) # dividing each by its total genes
  countHitPathways <- countHitPathways[order(countHitPathways[, 4], decreasing = TRUE),] # resorting with decreasing normalized values to decrease execution time in the following monte carlo for loop
  # assign transition probability values to pathways according their gene participation values
  sum_hits <- sum(as.numeric(as.character(countHitPathways[, 4]))) # 4th column with normalized probability values
  if (sum_hits == 0){
    pI <- sample(1:nP, 1) # random pathway restart in the case of non participating genes
    # cat(sprintf("Gene(s) not found in any pathway. New random pathway: %d #\n", pI))
  } else{
    rand <- runif(1, min = 0, max = sum_hits) # using random uniform distribution due to decimal values
    pathway_selector <- 0
    for (i in 1:nP){
      pathway_selector <- pathway_selector + as.numeric(countHitPathways[i, 4])
      if (rand <= pathway_selector){
        pI <- as.numeric(V(pG)[countHitPathways[i, 1]]) # getting id of chosen pathway
        # cat(sprintf("Gene(s) exist in %d pathway(s). New random -by score- pathway: %d, %s\n", colSums(countHitPathways != 0), pI, V(pG)[pI]$name))
        break
      }
    }
  }
  # pathway chosen, begin to walk the shortest path from current to next pathway
  shortestPathString <- all_shortest_paths(pG, current_pI, pI) # find all possible shortest paths between these two 
  sh_p_iter <- 1 # default, in case of one shortest path
  if (length(shortestPathString[[1]]) > 1){ # if more than one shortest paths
    sh_p_iter <- sample(1:length(shortestPathString[[1]]), 1) # choose random among equal shortest paths
    # cat(sprintf("%d paths between %s and %s\n", length(shortestPathString[[1]]), V(pG)$name[current_pI], V(pG)$name[pI]))
  }
  strSplit <- strsplit(shortestPathString[[1]][[sh_p_iter]]$name, "\t")
  # firstly increment the rankedPathways scores
  # cat(print(as.character(strSplit)))
  rankedPathways[as.character(strSplit), ] <<- rankedPathways[as.character(strSplit), ] + 1 # global variable assignment
  # then increment the rankedEdges scores
  num_walked_pathways <- length(strSplit)
  # cat(sprintf("Length of chosen shortest path is: %d #\n", num_walked_pathways))
  if (num_walked_pathways > 1){ # if walker did not stay at same pathway
    for (i in 2:num_walked_pathways){ # for all walked edges (couples of pathways)
      first_search <- which(rankedEdges[, 1] %in% strSplit[[i-1]]) %in% which(rankedEdges[, 2] %in% strSplit[[i]])
      if (any(first_search)){
        rankedEdges[which(rankedEdges[, 1] %in% strSplit[[i-1]])[which(first_search)], 3] <<- as.double(rankedEdges[which(rankedEdges[, 1] %in% strSplit[[i-1]])[which(first_search)], 3]) + 1
        # cat(print(rankedEdges[which(rankedEdges[, 1] %in% strSplit[[i-1]])[which(first_search)], ]))
      } else{
        second_search <- which(rankedEdges[, 2] %in% strSplit[[i-1]]) %in% which(rankedEdges[, 1] %in% strSplit[[i]])
        if(any(second_search)){
          rankedEdges[which(rankedEdges[, 2] %in% strSplit[[i-1]])[which(second_search)], 3] <<- as.double(rankedEdges[which(rankedEdges[, 2] %in% strSplit[[i-1]])[which(second_search)], 3]) + 1
          # cat(print(rankedEdges[which(rankedEdges[, 2] %in% strSplit[[i-1]])[which(second_search)], ]))
        }
      }
    }
  }
  return(pI)
}

# check similarity of two matrices (used in order to decide convergence)
check_matrix_similarity <- function(a, b){
  similarity_metric <- colSums(a == b)/length(a) # length(a) == length(b)
  return(similarity_metric) # [0, 1]
}

# convergence function based on n steps or stability of results based on a converge factor if n == 0
converged <- function(c){
  converge <-  FALSE
  rankedPathways <<- as.matrix(rankedPathways[order(rankedPathways, decreasing = TRUE),1]) # sorting ranked pathways decreasingly
  similarity <- check_matrix_similarity(as.matrix(rownames(rankedPathways)), as.matrix(rownames(previousRankedPathways))) # calculate similarity of current and last ranked-position results
  lastSimilarityIndexes <<- rbind(matrix(lastSimilarityIndexes[2:number_of_last_variances_to_check,]), similarity) # keep only last number_of_last_variances_to_check
  converge_variance <- var(lastSimilarityIndexes)
  cat(sprintf("Step: %d | Similarity: %f | Variance: %f | ", c, similarity, converge_variance))
  # writing outputs
  write(converge_variance, "converge_variance.txt", append = TRUE)
  write(similarity, "converge_results.txt", append = TRUE)
  end_time <<- Sys.time()
  dif <- end_time - start_time
  cat(sprintf("Time: %f\n", dif))
  saveRDS(rankedPathways, "rankedPathways.rds")
  saveRDS(rankedEdges, "rankedEdges.rds")
  start_time <<- Sys.time()
  if (similarity >= converging_factor && converge_variance <= converging_variance) converge <- TRUE
  previousRankedPathways <<- rankedPathways
  return(converge)
}

infOmicsWalk <- function(gGraph, pGraph) {
  num_genes <- length(V(gGraph))
  num_pathways <- length(V(pGraph))
  gene_iterator <- sample(1:num_genes, 1) # assign starting gene to gene_iterator
  pathway_iterator <- sample(1:num_pathways, 1) # assign starting pathway to pathway_iterator
  cat(sprintf("~ PathWalks with %d map genes and %d pathways.\n~ Starting at gene %s and at pathway %s.\n", num_genes, num_pathways, V(gGraph)[gene_iterator]$name, V(pGraph)[pathway_iterator]$name))
  iteration <- 0
  start_time <<- Sys.time()
  repeat{ # until converged
    iteration <- iteration + 1
    # cat(sprintf("%d\n", iteration))
    # gene map steps
    geneStepResults <- geneStepMultipleRWR(gGraph, gene_iterator, iteration) # iteration for restart
    gene_iterator <- geneStepResults$gene_it
    # pathway network steps
    pathway_iterator <- pathwayStepMultiple(pGraph, geneStepResults$gene, pathway_iterator)
    if ((iteration %% check_similarity_timer) == 0)if (converged(iteration)) break # iteration div timer alarms to check for convergence
  }
  return(1)
}

randGenerator <- function() { #generates a random alphanumeric string
  a <- do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE)) #1 string with 5 characters
  paste0(a, sprintf("%04d", sample(9999, 1, TRUE)), sample(LETTERS, 1, TRUE)) #paste with 4digits and 1 last character
}

translateNames <- function(m){
  for (i in 1:length(m)){ # for every pathway
    name <- m[i]
    m[i] <- pathway_translation_matrix[which(pathway_translation_matrix[, 1] %in% name), 2]
  }
  return(m)
}

printOutputRankedPathways <- function(rand){
  translatedMatrix <- translateNames(rownames(rankedPathways))
  translated_rankedPathways <- cbind(translatedMatrix, rankedPathways[, 1])
  outName <- paste(rand, "rankedPathways.tsv", sep = "_")
  unlink(outName) # if file already exists
  write.table(translated_rankedPathways, outName, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  return(1)
}

printOutputEdgelist <- function(rand){
  translatedMatrix1 <- translateNames(rankedEdges[, 1])
  translatedMatrix2 <- translateNames(rankedEdges[, 2])
  translated_rankedEdges <- cbind(translatedMatrix1, translatedMatrix2, rankedEdges[, 3])
  translated_rankedEdges <- as.data.frame(translated_rankedEdges)
  translated_rankedEdges[, 3] <- as.numeric(as.character(translated_rankedEdges[, 3])) # coerce the relevant column to numeric
  translated_rankedEdges <- translated_rankedEdges[order(as.numeric(translated_rankedEdges[, 3]), decreasing = TRUE), ]
  outName <- paste(rand, "rankedEdges.tsv", sep = "_")
  unlink(outName) # if file already exists
  write.table(translated_rankedEdges, outName, row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
  return(1)
}

louvainClustering <- function(filename){
  edgelist <- as.matrix(read.delim(filename, header = FALSE))
  graph <- graph_from_edgelist(edgelist[, 1:2], directed = FALSE)
  E(graph)$weight <- as.double(edgelist[, 3])
  # remove loops and multiple edges, simplify sum aggregates same edges
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "sum"))
  louvainClusters <- cluster_louvain(graph)
  outFile <- paste(randStr, "clusters.txt", sep="_")
  if (file.exists(outFile)) file.remove(outFile)
  for (i in 1:length(louvainClusters)){
    write.table(paste("\nCluster ", i, " - ", length(louvainClusters[[i]]), " pathways:\n", sep=""), outFile, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
    write.table(louvainClusters[i], outFile, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  }
  return(louvainClusters)
}

#####

# MAIN ####
# args
args = commandArgs(trailingOnly = TRUE) #allow use of args
if (length(args) == 0) {
  stop("At least two file arguments (gene and pathway networks) must be supplied.", call.=FALSE)
} else if (length(args) == 1) {
  stop("At least two file arguments (gene and pathway networks) must be supplied.", call.=FALSE)
} else if (length(args) == 2){
  # defaults
  gene_network_file <- args[1] # "geneEdgelistIPF.tsv" # "geneEdgelistAD.tsv"
  pathway_network_file <- args[2] # "hsaPathwayEdgelistIPF.tsv" # "hsaPathwayEdgelistAD.tsv" # "hsaPathwayEdgelist_commonGenes.tsv"
  converging_factor <- 0.95 # similarity threshold between last results needed to converge, values greater than 1 will allow the algorithm to run infinitely, higher values greatly increase execution time
  gene_restart_timer <- 50 # iterations until walker restarts at a random node of the gene map
  check_similarity_timer <- 100 # iterations until the converging criteria are checked by the algorithm
  number_of_last_variances_to_check <- 10 # variable that helps avoiding algorithm quits due to randomly similar results
  converging_variance <- 0.003 # variable that helps avoiding algorithm quits due to randomly similar results 
} else if (length(args) == 3){
  gene_network_file <- args[1]
  pathway_network_file <- args[2]
  converging_factor <- as.numeric(args[3])
  gene_restart_timer <- 50
  check_similarity_timer <- 100
  number_of_last_variances_to_check <- 10
  converging_variance <- 0.003
} else if (length(args) == 4){
  gene_network_file <- args[1]
  pathway_network_file <- args[2]
  converging_factor <- as.numeric(args[3])
  gene_restart_timer <- as.numeric(args[4])
  check_similarity_timer <- 100
  number_of_last_variances_to_check <- 10
  converging_variance <- 0.003
} else if (length(args) == 5){
  gene_network_file <- args[1]
  pathway_network_file <- args[2]
  converging_factor <- as.numeric(args[3])
  gene_restart_timer <- as.numeric(args[4])
  check_similarity_timer <- as.numeric(args[5])
  number_of_last_variances_to_check <- 10
  converging_variance <- 0.003
} else if (length(args) == 6){
  gene_network_file <- args[1]
  pathway_network_file <- args[2]
  converging_factor <- as.numeric(args[3])
  gene_restart_timer <- as.numeric(args[4])
  check_similarity_timer <- as.numeric(args[5])
  number_of_last_variances_to_check <- as.numeric(args[6])
  converging_variance <- 0.003
} else if (length(args) == 7){
  gene_network_file <- args[1]
  pathway_network_file <- args[2]
  converging_factor <- as.numeric(args[3])
  gene_restart_timer <- as.numeric(args[4])
  check_similarity_timer <- as.numeric(args[5])
  number_of_last_variances_to_check <- as.numeric(args[6])
  converging_variance <- as.numeric(args[7])
}

# load library
cat(sprintf("Loading igraph library..\n"))
suppressMessages(library("igraph"))

# load translation files
cat(sprintf("Loading KEGG gene and pathway translation and linkage  files..\n"))
gene_translation_matrix <- as.matrix(read.delim("gene_translation_dec_2019.tsv", header = FALSE))
pathway_translation_matrix <- as.matrix(read.delim("hsa_pathway_tranlsation_dec_2019.tsv", header = FALSE))
gene_pathway_links_matrix <- as.matrix(read.delim("gene_pathway_links_dec_2019.tsv", header = FALSE))
# parse for pathways selection normalization
num_genes_per_pathway <- as.matrix(table(gene_pathway_links_matrix[,1]))
num_genes_per_pathway <- cbind(rownames(num_genes_per_pathway), num_genes_per_pathway)
colnames(num_genes_per_pathway) <- c("pathwayID", "total_genes")
rownames(num_genes_per_pathway) <- c()

# generate gene map
cat(sprintf("Generating Gene Map..\n")) # gene_network_file is args[1]
hsaGeneEdgelist <- translateInputGenes(gene_network_file)
geneGraph <- createGeneNetwork(hsaGeneEdgelist)

# generate pathway network
cat(sprintf("Generating Pathway Network..\n")) # pathway_network_file is args[2]
hsaPathwayEdgelist <- as.matrix(read.delim(pathway_network_file, header = FALSE))
pathwaysGraph <- createPathwayNetwork(hsaPathwayEdgelist)
# check if graph is fully connected, otherwise remove subgraphs
if (!is_connected(pathwaysGraph)) pathwaysGraph <- removeDetachedNodes(pathwaysGraph)

# Initialization of execution variables 
cat(sprintf("Initializing execution variables..\n"))
# initilization of matrix (and its previous version) with ordered pathway results
# the rankedPathways matrix will keep the scores for every pathway visited after walking an edge or staying on the same pathway
pathway_nodes <- V(pathwaysGraph)$name
rankedPathways <- matrix(0, nrow = length(pathway_nodes))
rownames(rankedPathways) <- pathway_nodes
previousRankedPathways <- rankedPathways
# removing non-participating pathways from num_genes_per_pathway, to replace a later merge with cbind for time purposes
pathway_nodes_matrix <- as.matrix(pathway_nodes)
pathway_nodes_matrix <- as.matrix(pathway_nodes_matrix[order(pathway_nodes_matrix[,1]),])
colnames(pathway_nodes_matrix) <- "pathwayID"
num_genes_per_pathway <- merge(pathway_nodes_matrix, num_genes_per_pathway)

# initilization of edge score output matrix
rankedEdges <- hsaPathwayEdgelist[, 1:2]
zero_column <- matrix(0, nrow = nrow(rankedEdges), ncol = 1)
rankedEdges <- cbind(rankedEdges, zero_column)

# initialization of the similarities matrix to check variance of values
lastSimilarityIndexes <- matrix(-1, nrow = number_of_last_variances_to_check , ncol = 1)

# initialization of global timer variables
start_time <- 0
end_time <- 0

# Execution start
cat(sprintf("Starting the algorithm execution..\n"))
# # read from file, if previously stopped mid-execution
# rankedPathways <- readRDS("rankedPathways.rds")
# rankedEdges <- readRDS("rankedEdges.rds")
start <- Sys.time()
done <- suppressWarnings(infOmicsWalk(geneGraph, pathwaysGraph))
end <- Sys.time()
dif <- end - start
cat(sprintf("Total PathWalks execution time: %f\n", dif))

# writing outputs; pathway ranks and pathway edgelist with scores
cat(sprintf("Printing Outputs..\n"))
randStr <- randGenerator() #random string for unique name
done <- printOutputRankedPathways(randStr)
done <- printOutputEdgelist(randStr)

# Louvain community detection
cat(sprintf("Printing Communities as detected by the Louvain method..\n"))
louvainClusters <- louvainClustering(paste(randStr, "rankedEdges.tsv", sep="_"))

#####
