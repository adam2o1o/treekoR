#' runHOPACH
#'
#' @param data dataframe containing the median expression of the clusters/cell
#' types
#' @param K  positive integer specifying the maximum number of levels in the
#' tree. Must be 15 or less, due to computational limitations (overflow)
#' @param kmax integer between 1 and 9 specifying the maximum number of children
#' at each node in the tree
#' @param dissimilarity_metric metric used to calculate dissimilarities
#' between clusters/cell types
#'
#' @import hopach
#' @return a list containing the groups each cluster belongs to at each level of
#' the hopach tree
#'
#' @examples
#' library(SingleCellExperiment)
#' library(data.table)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_med_dt <- as.data.table(exprs)
#' clust_med_dt[, cluster_id := clusters]
#' res <- clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
#' res2 <- res[,.SD, .SDcols = !c('cluster_id')]
#'
#' hopach_res <- runHOPACH(as.data.frame(scale(res2)))
#' @export
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
#' @importFrom ape read.tree
#' @return a phylogram converted from the outputted list from the runHOPACH
#' function
#' @examples
#' library(SingleCellExperiment)
#' library(data.table)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_med_dt <- as.data.table(exprs)
#' clust_med_dt[, cluster_id := clusters]
#' res <- clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
#' res2 <- res[,.SD, .SDcols = !c('cluster_id')]
#'
#' hopach_res <- runHOPACH(as.data.frame(scale(res2)))
#' phylo <- hopachToPhylo(hopach_res)
#'
#' @export
hopachToPhylo <- function(res) {
  cutree_list_df <- do.call(cbind, res$cutree_list)
  for (i in seq_len(ncol(cutree_list_df)-1)) {
    cutree_list_df[,i] <- paste0(i,"_",cutree_list_df[,i])
  }

  # We add < & > to each cluster for string matching purposes
  base_nodes <- paste0("<",as.vector(cutree_list_df[,ncol(cutree_list_df)]),">")
  # A vector containing the current grouping of nodes as the loop iterates
  current_nodes <- base_nodes

  # This loop is to iterate through the levels of the tree and
  # sequentially combine clusters in a format to be read into a
  # phylogenetic tree
  # ------------------------------------------------------------
  for (i in rev(seq_len(ncol(cutree_list_df)-1))) {
    cur_lev_nodes <- as.vector(cutree_list_df[,i])
    cur_lev_n_t <- table(cur_lev_nodes)
    # Get parent nodes with more than two children
    unique_nodes <- names(cur_lev_n_t[which(cur_lev_n_t >= 2)])
    # For each parent node, we combine the corresponding children nodes
    for (node in unique_nodes) {
      base_nodes_combo <- base_nodes[which(cur_lev_nodes== node)]
      base_nodes_rgx_s <- paste(base_nodes_combo, collapse="|")
      if (sum(grepl(base_nodes_rgx_s, current_nodes)) > 1) {
        nodes_to_combine <- current_nodes[grepl(base_nodes_rgx_s, current_nodes)]

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


#' getClusterTree
#' This function takes a CATALYST sce with clusters and creates a hierarchical tree
#'
#' @param exprs a dataframe containing single cell expression data
#' @param clusters a vector representing the cell type or cluster of each cell
#' (can be character or numeric). If numeric, cluster names need to be consecutive
#' starting from 1.
#' @param hierarchy_method a string indicating the hierarchical tree construction
#'  method to be used
#' @param hopach_kmax integer between 1 and 9 specifying the maximum number of
#' children at each node in the tree
#' @param hopach_K positive integer specifying the maximum number of levels in the
#' tree. Must be 15 or less, due to computational limitations (overflow)
#' @param scale_exprs boolean indicating whether to scale median cluster expression
#' data before constructing hierarchical tree
#'
#' @return
#' @import data.table
#' @importFrom ggtree ggtree
#' @importFrom ape as.phylo
#' @importFrom stats dist hclust median
#'
#' @return a list containing the cluster median frequencies and a phylogram of the
#' hierarchical tree
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
getClusterTree <- function(exprs,
                           clusters,
                           hierarchy_method="hopach",
                           hopach_kmax = 5,
                           hopach_K = 10,
                           scale_exprs=TRUE) {
  clust_med_dt <- as.data.table(exprs)
  clust_med_dt[, cluster_id := clusters]
  # data table containing median
  res <- clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
  res2 <- res[,.SD, .SDcols = !c('cluster_id')]
  rownames(res2) <- res[["cluster_id"]]

  if (scale_exprs) {
    res_unscaled <- res2
    res2[, (colnames(res2)) := lapply(.SD, scale), .SDcols=colnames(res2)]
  } else {
    res_unscaled <- res2
  }

  if (hierarchy_method == "hopach") {
    hp_dend <- runHOPACH(data = as.data.frame(res2),
                         kmax=hopach_kmax,
                         K=hopach_K)
    hc_phylo <- hopachToPhylo(hp_dend)
    hc_phylo$tip.label <- rownames(res2)[as.numeric(hc_phylo$tip.label)]
  } else {
    clust_dist <- dist(res2)
    hc_dend <- hclust(clust_dist, method=hierarchy_method)
    hc_phylo <- as.phylo(hc_dend)
    hc_phylo$tip.label <- as.character(res[["cluster_id"]])
  }

  return(list(
    median_freq = res_unscaled,
    clust_tree = hc_phylo
  ))
}

#' Title
#'
#' @param tree a ggtree object
#'
#' @return a ggtree object with the data containing a column with the clusters
#' contained in each node
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
  d$clusters <- lapply(d$clusters, unlist)
  d$clusters <- lapply(d$clusters, unique)
  d$clusters <- lapply(d$clusters, function(x) x[!is.na(x)])

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
#' @param paired a boolean indicating whether to performed paired t-tests
#' (not yet tested)
#' @param phylo a ggtree object
#' @param clusters a vector representing the cell type or cluster of each cell
#' (can be character or numeric). If numeric, cluster names need to be consecutive
#' starting from 1.
#' @param classes a vector containing the patient outcome/class each cell belongs to
#' @param samples a vector identifying the patient each cell belongs to
#' @param pos_class_name a character indicating which class is positive
#' @param subjects a vector containing which subject the cell belongs to, used
#' to identify matched samples in paired t-tests (not yet tested)
#'
#' @importFrom stats t.test
#'
#' @return a ggtree object with significance testing results in embedded data
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
#'
#' tested_tree <- testTree(clust_tree$clust_tree,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes,
#'                         pos_class_name=NULL,
#'                         subjects=NULL,
#'                         paired = FALSE)
testTree <- function(phylo,
                     clusters,
                     samples,
                     classes,
                     pos_class_name=NULL,
                     subjects=NULL,
                     paired = FALSE){
  if (!is.null(pos_class_name)) {
    if (!pos_class_name %in% unique(classes)) {
      stop("'pos_class_name' needs to be in classes")
    }
  }
  if (length(unique(classes)) > 2) {
    stop("length(unique(classes)) > 2.
         treekoR can currently only test between two classes.")
  }
  t <- findChildren(ggtree(phylo, branch.length="none"))
  td <- t$data
  if(paired == TRUE){
    samp2Group <- unique(data.frame(subjects, samples, classes))
    samp2Group <- samp2Group[samp2Group$subjects %in%
                               names(which(table(samp2Group$subjects)==2)),]
    ### Put catch error in here
    groupA <- as.character(samp2Group[samp2Group$classes==pos_class_name,'samples'])
    groupB <- as.character(samp2Group[samp2Group$classes==neg_class_name,'samples'])
  } else {
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
  for( i in seq_len(nrow(td))){
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
      testAll <- t.test(ratioAll[groupA],ratioAll[groupB],paired = paired)
      statAll[i] <- testAll$statistic
      pvalAll[i] <- testAll$p.value
      testParent <- t.test(ratioParent[groupA],ratioParent[groupB],paired = paired)
      statParent[i] <- testParent$statistic
      pvalParent[i] <- testParent$p.value
    }
  }
  td[, c("statAll", "statParent", "pvalAll", "pvalParent")] <-
    data.frame(statAll, statParent, pvalAll, pvalParent)
  t$data <- td
  return(t)
}

#' getCellProp
#'
#' @param phylo a phylogram with tip.labels corresponding to cell types/cluster
#' contained in 'clusters' vector
#' @param clusters a vector representing the cell type or cluster of each
#' cell (can be character or numeric). If numeric, cluster names need to be
#' consecutive starting from 1.
#' @param samples a vector identifying the patient each cell belongs to
#' @param classes a vector containing the patient outcome/class each cell belongs to
#'
#' @return a dataframe containing proportions calculated for each sample
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
#'
#' prop_df <- getCellProp(clust_tree$clust_tree,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes)
getCellProp <- function(phylo,
                        clusters,
                        samples,
                        classes){
  t <- findChildren(ggtree(phylo, branch.length="none"))
  td <- t$data
  td$parentClusters <- td$clusters[match(td$parent, td$node)]
  # Remove root node
  td <- td[-which(td$node == td$parent),]

  # Get dataframe of parent proportions
  prop_par <- data.frame(mapply(function(x,y) {
    tapply(clusters %in% unlist(x),
           samples, sum) /tapply(clusters %in% unlist(y), samples, sum)
  }, td$clusters, td$parentClusters))

  colnames(prop_par) <- paste0("prop_parent_", ifelse(is.na(td$label),
                                                      td$node, td$label))

  # Get dataframe of absolute proportions
  prop_all <- data.frame(mapply(function(x) {
    tapply(clusters %in% unlist(x), samples, sum)/table(samples)
  }, td$clusters))

  colnames(prop_all) <- paste0("prop_all_", ifelse(is.na(td$label),
                                                   td$node, td$label))

  # Mapping from sample to class outcome
  samp_class <- unique(cbind(as.character(samples), as.character(classes)))

  prop_df <- cbind(data.frame(sample = names(table(samples)),
                              class = samp_class[,2][match(names(table(samples)),
                                                           samp_class[,1])]),
                   prop_all,
                   prop_par)

  return(prop_df)
}

#' getTreeResults
#'
#' @param testedTree a ggtree object outputed from testTree()
#' @param sort_by whether to sort by p-values testing via proportions to parent or
#' p-values testing via absolute proportions. Values can can be c(NA, "parent", "all")
#'
#' @return a dataframe with hierarchical tree nodes, corresponding clusters and
#' corresponding significance testing results
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(COVIDSampleData)
#'
#' sce <- DeBiasi_COVID_CD8_samp
#' exprs <- t(assay(sce, "exprs"))
#' clusters <- colData(sce)$cluster_id
#' classes <- colData(sce)$condition
#' samples <- colData(sce)$sample_id
#'
#' clust_tree <- getClusterTree(exprs,
#'                              clusters,
#'                              hierarchy_method="hopach")
#'
#' tested_tree <- testTree(clust_tree$clust_tree,
#'                         clusters=clusters,
#'                         samples=samples,
#'                         classes=classes,
#'                         pos_class_name=NULL,
#'                         subjects=NULL,
#'                         paired = FALSE)
#'
#' res_df <- getTreeResults(tested_tree)
#'
#' head(res_df, 10)
getTreeResults <- function(testedTree,
                           sort_by = "parent"){
  if (!sort_by %in% c("parent", "all")) {
    stop("'sort_by' must be 'parent' or 'all'")
    }

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
