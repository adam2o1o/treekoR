#' Title
#'
#' @param data dataframe containing the median expression of the clusters/cell types
#' @param K  positive integer specifying the maximum number of levels in the tree. Must be
#' 15 or less, due to computational limitations (overflow)
#' @param kmax integer between 1 and 9 specifying the maximum number of children at each
#' node in the tree
#' @param dissimilarity_metric metric used to calculate dissimilarities between clusters/cell types
#'
#' @return
#' @export
#' @import hopach
runHOPACH <- function(data, K = 10, kmax = 5, dissimilarity_metric = "cor") {
  # Compute the distance matrix using correlation
  dist <- distancematrix(data, d = dissimilarity_metric)
  # Run HOPACH function
  clustresult <- hopach(data, K = K , dmat = dist, kmax = kmax)
  # Rearrange HOPACH results
  final_labels <- strsplit(as.character(format(clustresult$final$labels,
                                               scientific = FALSE)), "")
  # Get each level seperate
  level_list <- list()
  for (i in seq_len(nchar(clustresult$final$labels[1]))) {
    if (i != 1) {
      level_list[[i]] <- paste(level_list[[i - 1]],
                               unlist(lapply(final_labels, "[[", i)), sep = "")
    }else{
      level_list[[i]] <- unlist(lapply(final_labels, "[[", i))
    }
  }
  # Generate cutree results
  cutree_list <- list()
  # First level (All ones)
  cutree_list[[1]] <- rep(1, length(rownames(data)))
  names(cutree_list[[1]]) <- rownames(data)
  cutree_list[2:(length(level_list) + 1)] <- lapply(level_list, function(x) {
    x <- as.numeric(as.factor(x))
    names(x) <- rownames(data)
    x
  })
  # Last levels will be all distinct numbers
  cutree_list[[length(cutree_list) + 1]] <- seq_len(length(clustresult$final$labels))
  names(cutree_list[[length(cutree_list)]]) <- rownames(data)

  return(list(cutree_list = cutree_list))
}

#' Title
#'
#' @param res an object returned from the runHOPACH() function
#'
#' @return
#' @importFrom ape read.tree
#'
#' @export
hopachToPhylo <- function(res) {

  cutree_list_df <- do.call(cbind, res$cutree_list)
  for (i in 1:(ncol(cutree_list_df)-1)) {
    cutree_list_df[,i] <- paste0(i,"_",cutree_list_df[,i])
  }

  # We add < & > to each cluster for string matching purposes
  base_nodes <- paste0("<",as.vector(cutree_list_df[,ncol(cutree_list_df)]),">")
  # A vector containing the current grouping of nodes as the loop iterates
  current_nodes <- base_nodes

  #' This loop is to iterate through the levels of the tree and
  #' sequentially combine clusters in a format to be read into a
  #' phylogenetic tree
  #' ------------------------------------------------------------
  for (i in (ncol(cutree_list_df)-1):1) {
    # Get parent nodes with more than two children
    unique_nodes <- names(table(cutree_list_df[,i])[which(table(cutree_list_df[,i]) >= 2)])
    # For each parent node, we combine the corresponding children nodes
    for (node in unique_nodes) {
      base_nodes_combo <- base_nodes[which(as.vector(cutree_list_df[,i])== node)]

      if (sum(grepl(paste(base_nodes_combo, collapse="|"), current_nodes)) > 1) {
        nodes_to_combine <- current_nodes[grepl(paste(base_nodes_combo, collapse="|"), current_nodes)]

        current_nodes <- c(current_nodes[!current_nodes %in% nodes_to_combine],
                           paste0("(", paste(nodes_to_combine, collapse=","), ")")
        )
      }
    }
  }

  tree_string <- paste0(gsub("<|>", "", current_nodes), ";")
  # Conver tree string into a phylo object
  tree <- read.tree(text = tree_string)

  return(tree)
}


#' Title
#' This function takes a CATALYST sce with clusters and creates a hierarchical tree
#'
#' @param exprs a dataframe containing single cell expression data
#' @param clusters a vector representing the cell type or cluster of each cell (can be character or numeric)
#' @param hierarchy_method a string indicating the hierarchical tree construction method to be used
#' @param hopach_kmax integer between 1 and 9 specifying the maximum number of children at each
#' node in the tree
#' @param hopach_K positive integer specifying the maximum number of levels in the tree. Must be
#' 15 or less, due to computational limitations (overflow)
#'
#' @return
#' @import data.table
#' @importFrom ggtree ggtree
#' @importFrom ape as.phylo
#'
#' @export
getClusterTree <- function(exprs,
                           clusters,
                           hierarchy_method="hopach",
                           hopach_kmax = 5,
                           hopach_K = 10) {
  clust_med_dt <- as.data.table(exprs)
  clust_med_dt[, cluster_id := clusters]
  # data table containing median
  res <-clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
  res2 <- res[,.SD, .SDcols = !c('cluster_id')]
  rownames(res2) <- res[["cluster_id"]]

  if (hierarchy_method == "hopach") {
    hp_dend <- runHOPACH(data = as.data.frame(scale(res2)),
                         kmax=hopach_kmax,
                         K=hopach_K)
    hc_phylo <- hopachToPhylo(hp_dend)
    hc_phylo$tip.label <- rownames(res2)[as.numeric(hc_phylo$tip.label)]
  } else {
    clust_dist <- dist(scale(res2))
    hc_dend <- hclust(clust_dist, method=hierarchy_method)
    hc_phylo <- as.phylo(hc_dend)
    hc_phylo$tip.label <- as.character(res[["cluster_id"]])
  }

  return(list(
    median_freq = res2,
    clust_tree = hc_phylo
  ))
}

#' Title
#'
#' @param tree a ggtree object
#'
#' @return
findChildren <- function(tree) {
  d <- tree$data
  d$clusters <- d$label
  d$clusters <- as.list(as.character(d$clusters))
  uNodes <- sort(unique(d$x), decreasing = TRUE)
  for(x in uNodes){
    nodes = as.matrix(d[which(d$x==x), "node"])
    for(n in nodes){
      parentNode <- as.character(d$parent[which(d$node == n)])
      childClusters <- d$clusters[[which(d$node == n)]]
      parentNodeInd <- which(d$node == parentNode)
      if (is.na(d$clusters[parentNodeInd])) {
        d$clusters[[parentNodeInd]] <- childClusters
      } else {
        d$clusters[[parentNodeInd]] <- c(d$clusters[[parentNodeInd]], childClusters)
      }
    }
  }
  d$clusters <- sapply(d$clusters, unlist)
  d$clusters <- lapply(d$clusters , unique)
  d$clusters <- lapply(d$clusters , function(x) x[!is.na(x)])

  tree$data <- d
  return(tree)
}

#' Title
#' @description This function takes a hierarchical tree of the cluster
#' medians of a cytometry dataset, and then uses this structure to perform
#' t-tests between conditions of patients testing for difference using the
#' proportion of cluster relative to sample's n and proportion of cluster
#' relative to sample's n of hierarchical parent cluster.
#' Takes a ggtree object and returns a ggtree object with testing results
#' appended in the data
#'
#' @param paired a boolean indicating whether to performed paired t-tests (not yet tested)
#' @param phylo a ggtree object
#' @param exprs a dataframe containing the clusters for each cell, the
#' sample id, the subject id (needed for paired tests), and the group
#' that the subject/sample belongs to
#' @param clusters a vector representing the cell type or cluster of each cell (can be character or numeric)
#' @param classes a vector containing the patient outcome/class each cell belongs to
#' @param samples a vector identifying the patient each cell belongs to
#' @param pos_class_name a character indicating which class is positive
#' @param subjects a vector containing which subject the cell belongs to, used
#' to identify matched samples in paired t-tests (not yet tested)
#'
#' @return
#' @export
testTree <- function(phylo,
                     exprs,
                     clusters,
                     samples,
                     classes,
                     pos_class_name=NULL,
                     subjects=NULL,
                     paired = FALSE){
  t <- ggtree(phylo, branch.length="none")
  t <- findChildren(t)
  td <- t$data

  if(paired == TRUE){
    samp2Group <- unique(data.frame(subjects, samples, classes))
    samp2Group <- samp2Group[samp2Group$subjects%in%names(which(table(samp2Group$subjects)==2)),]
    ### Put catch error in here
    groupA <- as.character(samp2Group[samp2Group$classes==pos_class_name,'samples'])
    groupB <- as.character(samp2Group[samp2Group$classes==neg_class_name,'samples'])

  }else{
    if (is.null(pos_class_name)) {
      pos_class_name <- unique(classes)[1]
      neg_class_name <- unique(classes)[2]
    } else {
      neg_class_name <- unique(classes)[unique(classes) != pos_class_name]
    }
    samp2Group <- unique(data.frame(samples, classes))
    groupA <- as.character(samp2Group[samp2Group$classes==pos_class_name,'samples'])
    groupB <- as.character(samp2Group[samp2Group$classes==neg_class_name,'samples'])
  }
  statParent <- statAll <- pvalParent <- pvalAll <- NULL

  for( i in 1:nrow(td)){
    childClusters <- td[i,'clusters'][[1]][[1]]
    parent <- td[i,'parent'][[1]]
    parentClusters <- td[td$node == parent,'clusters'][[1]][[1]]

    if(mean(parentClusters %in% childClusters)==1){
      statAll[i] <- 0
      pvalAll[i] <- 1
      statParent[i] <- 0
      pvalParent[i] <- 1
    } else {
      numAll <- table(samples)
      numParent <- tapply(clusters %in% parentClusters, samples , sum)
      numNode <- tapply(clusters %in% childClusters, samples, sum)

      ratioAll <- numNode/pmax(numAll[names(numNode)],1)
      ratioParent <- numNode/pmax(numParent[names(numNode)],1)

      if(paired == TRUE){
        testAll <- t.test(ratioAll[groupA],ratioAll[groupB],paired = TRUE)
        statAll[i] <- testAll$statistic
        pvalAll[i] <- testAll$p.value


        testParent <- t.test(ratioParent[groupA],ratioParent[groupB],paired = TRUE)
        statParent[i] <- testParent$statistic
        pvalParent[i] <- testParent$p.value
      } else {
        testAll <- t.test(ratioAll[groupA],ratioAll[groupB])
        statAll[i] <- testAll$statistic
        pvalAll[i] <- testAll$p.value


        testParent <- t.test(ratioParent[groupA],ratioParent[groupB])
        statParent[i] <- testParent$statistic
        pvalParent[i] <- testParent$p.value
      }

    }
  }

  td$statAll <- statAll
  td$statParent <- statParent
  td$pvalAll <- pvalAll
  td$pvalParent <- pvalParent

  t$data <- td
  return(t)

}


#' Title
#'
#' @param testedTree a ggtree object outputed from testTree()
#' @param sort_by whether to sort by p-values testing via proportions to parent or
#' p-values testing via absolute proportions. Values can can be c(NA, "parent", "all")
#'
#' @return
#' @export
getTreeResults <- function(testedTree,
                           sort_by = "parent"){
  res <- as.data.frame(testedTree$data[,c(
    "parent", "node", "isTip",
    "clusters", "statAll", "statParent",
    "pvalAll", "pvalParent"
    )])

  if (sort_by == "parent") {
    res <- res[order(res$pvalParent),]
  }

  if (sort_by == "all") {
    res <- res[order(res$pvalAll),]
  }

  return(res)
}
